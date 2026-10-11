import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


def test_select_controls_matches_magnitude_and_sky_and_excludes_hosts():
    from calibrate_null_real import select_controls
    hosts = pd.DataFrame(dict(source_id=[1, 2], ra=[10.0, 200.0], dec=[0.0, -30.0],
                              phot_g_mean_mag=[11.0, 13.0]))
    rng = np.random.default_rng(0)
    n = 400
    cand = pd.DataFrame(dict(source_id=np.arange(100, 100 + n),
                             ra=np.r_[rng.uniform(8, 12, n // 2), rng.uniform(198, 202, n // 2)],
                             dec=np.r_[rng.uniform(-2, 2, n // 2), rng.uniform(-32, -28, n // 2)],
                             phot_g_mean_mag=rng.uniform(9, 15, n), ruwe=rng.uniform(0.8, 1.5, n)))
    cand.loc[0, "source_id"] = 1                                   # a host must never be its own control
    ctl = select_controls(cand, hosts, k=5, max_dg=0.5, max_sep_deg=3.0, max_ruwe=1.1, seed=1)
    assert set(ctl.host_source_id) == {1, 2}
    assert (ctl.groupby("host_source_id").size() <= 5).all()
    assert 1 not in set(ctl.source_id) and 2 not in set(ctl.source_id)
    assert (ctl.ruwe < 1.1).all()
    m = ctl.merge(hosts.rename(columns={"source_id": "host_source_id", "phot_g_mean_mag": "g_host",
                                        "ra": "ra_h", "dec": "dec_h"}), on="host_source_id")
    assert (np.abs(m.phot_g_mean_mag - m.g_host) <= 0.5).all()
    assert ctl.source_id.is_unique
