"""Independent cross-check: fit Gaia-4 (DR4 pre-release) with exogaia (T. Stolker, v0.1.3).

Run in a SEPARATE venv (exogaia pins healpy, which does not build on Windows):
    python -m venv exo-venv
    exo-venv/Scripts/python -m pip install --no-deps exogaia==0.1.3
    exo-venv/Scripts/python -m pip install astroquery beartype corner dynesty emcee h5py \
        kepler.py matplotlib nsstools numpy pandas pooch scipy tqdm
    # the 0.1.3 wheel omits its data tables; copy them from the v0.1.3 tag into
    # <exo-venv>/Lib/site-packages/data/: table_u0_g_c_p5.txt, holl2023_sigma.csv
    # put GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml in <HOME>/.exogaia/
    #   (sha256 f81f4dc1...6f3ce3, identical to ESA's file in data/raw/dr4_prerelease/)
    exo-venv/Scripts/python -I scripts/crosscheck_exogaia.py

healpy and imf are stubbed: exogaia only uses them for simulating data.
Result (2026-10-11), exogaia orbit_fit(inc_jitter=True): P = 573.07 +- 3.16 d, e = 0.292 +- 0.060,
a0 = 0.247 +- 0.008 mas, i = 121.77 +- 2.67 deg, omega = 180.6 +- 13.0 deg, Omega = 173.68 +- 2.59 deg,
jitter 0.049 mas, M_c = 1.010e-2 M_sun (14 function evaluations). Compare results/gaia4_validation.csv
and results/gaia4_period_profile.csv (our likelihood peaks at 579 d; 573 d is at -2 dlnL = 2.9).
"""
import sys
import types

for name in ("healpy", "imf", "imf.imf"):
    sys.modules[name] = types.ModuleType(name)
sys.modules["imf.imf"].make_cluster = None

from exogaia.data import GaiaAstrometry  # noqa: E402
from exogaia.leastsq import LeastSquares  # noqa: E402

ga = GaiaAstrometry(primary_mass=(0.644, 0.02), gaia_release="DR4", verbose=True)
ga.retrieve_dr4_prerelease(1457486023639239296, exclude_outliers=True, combine_ccds=False)
LeastSquares(ga).orbit_fit(inc_jitter=True, verbose=True)
