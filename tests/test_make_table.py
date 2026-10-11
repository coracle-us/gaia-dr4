import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


def _results():
    return pd.DataFrame(dict(
        hostname=["WASP-1", "HAT-P-2", "KELT-9"], source_id=[1, 2, 3],
        phot_g_mean_mag=[11.2, 8.6, 7.6], dist_pc=[380.0, 128.0, 206.0],
        classification=["orbit", "acceleration", "none"],
        dlnl2_orbit=[120.5, 60.0, 12.0], dlnl2_acc7=[3.0, 80.0, 1.0], dlnl2_acc9=[4.0, 95.0, 2.0],
        period_yr=[2.31, np.nan, np.nan], e=[0.21, np.nan, np.nan],
        post_mass_mj_p16=[7.1, np.nan, np.nan], post_mass_mj_p50=[8.0, np.nan, np.nan],
        post_mass_mj_p84=[9.2, np.nan, np.nan], post_inc_deg_p50=[88.0, np.nan, np.nan],
        post_min_mutual_inc_deg_p16=[0.5, np.nan, np.nan], post_min_mutual_inc_deg_p50=[2.0, np.nan, np.nan],
        post_min_mutual_inc_deg_p84=[5.0, np.nan, np.nan],
        m50_mj_P3=[12.0, 3.0, 9.0], m90_mj_P3=[19.0, 4.5, 14.0],
        flags=["twin_excluded_lit_rv;stability_confirmed", "twin_excluded_lit_rv", ""],
        status=["fitted"] * 3))


def test_latex_table_lists_detections_with_uncertainties():
    from make_table import latex_detection_table
    tex = latex_detection_table(_results())
    assert "\\begin{deluxetable" in tex and "\\end{deluxetable" in tex
    assert "WASP-1" in tex and "HAT-P-2" in tex and "KELT-9" not in tex
    assert "$8.0^{+1.2}_{-0.9}$" in tex                    # mass with asymmetric errors
    assert "acceleration" in tex.lower()


def test_machine_readable_table_contains_all_hosts(tmp_path):
    from astropy.table import Table
    from make_table import write_mrt
    path = tmp_path / "table.mrt"
    write_mrt(_results(), path)
    t = Table.read(path, format="ascii.mrt")
    assert len(t) == 3 and "m50_P3" in t.colnames and "Class" in t.colnames
