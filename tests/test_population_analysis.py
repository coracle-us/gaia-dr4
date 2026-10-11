import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


def _row(**kw):
    base = dict(hostname="H", source_id=1, classification="orbit", period_yr=2.0, post_mass_mj_p50=5.0,
                mass_mj=5.0, flags="twin_excluded_lit_rv;stability_confirmed")
    base.update(kw)
    return pd.Series(base)


def test_detection_selection():
    from population_analysis import is_detection
    assert is_detection(_row())
    assert not is_detection(_row(classification="acceleration"))
    assert not is_detection(_row(post_mass_mj_p50=40.0, flags="brown_dwarf_mass"))
    assert not is_detection(_row(period_yr=6.0))
    assert not is_detection(_row(flags="marginal"))


def test_population_script_runs(tmp_path):
    from population_analysis import main
    rng = np.random.default_rng(0)
    n = 60
    res = pd.DataFrame(dict(hostname=[f"H{i}" for i in range(n)], source_id=np.arange(n),
                            status="fitted", phot_g_mean_mag=rng.uniform(8, 13, n),
                            dist_pc=rng.uniform(30, 300, n), st_mass=1.0, s_ccd=0.02,
                            n_obs=800, pl_orbper=rng.uniform(1, 30, n), pl_orbincl=88.0,
                            classification="none", period_yr=np.nan, post_mass_mj_p50=np.nan,
                            mass_mj=np.nan, flags="", post_inc_deg_p50=np.nan,
                            post_inc_deg_p16=np.nan, post_inc_deg_p84=np.nan))
    for i in range(4):
        res.loc[i, ["classification", "period_yr", "post_mass_mj_p50", "post_inc_deg_p50",
                    "post_inc_deg_p16", "post_inc_deg_p84"]] = ["orbit", 2.0, 6.0, 88.0 + i, 86.0 + i, 90.0 + i]
    path = tmp_path / "results.csv"
    res.to_csv(path, index=False)
    main([str(path), "--out", str(tmp_path)])
    s = json.loads((tmp_path / "population_summary.json").read_text())
    occ = s["occurrence"]["all, aligned"]
    assert occ["n_detected"] == 4 and 0 < occ["median"] < 1
    assert "mutual_inclination_rayleigh_sigma_deg" in s
