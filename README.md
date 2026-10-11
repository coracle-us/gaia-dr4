# gaia-dr4-companions

Searching Gaia DR4 epoch astrometry (release: 2 Dec 2026) for outer giant companions to
stars hosting transiting giant planets. Target output: a Research Note of the AAS.

## Quick start (Windows, PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

See `CLAUDE.md` for the briefing, schedule, and rules.

## Pipeline

| Piece | Where |
|---|---|
| Epoch astrometry I/O (DR4 layout), DataLink / cache / simulated sources | `companions/epoch.py`, `dr4io.py`, `archive.py` |
| Single-star fit (reproduces ESA `gaiasupdate`), orbits, accelerations, transit-correlated noise | `companions/linear.py`, `orbit.py`, `accel.py`, `noise.py` |
| Orbit search, posteriors, joint astrometry + RV fit | `companions/search.py`, `posterior.py`, `joint.py` |
| DR4 sampling simulator, impostor physics, vetting | `companions/scanlaw.py`, `simulate.py`, `impostor.py`, `gaia_rv.py`, `vetting.py` |
| Per-host cascade (classification + vetting flags) | `companions/cascade.py` |
| Calibration (null, injection–recovery, stability, sensitivity) | `scripts/calibrate_null.py`, `injection_recovery.py`, `calibrate_stability.py`, `target_sensitivity.py`, `analyze_calibration.py` |
| Release day | `scripts/release_day.py`, `scripts/calibrate_null_real.py`, `scripts/make_table.py`, `config/release_day.json` |

Release-day procedure: [`docs/RELEASE_DAY.md`](docs/RELEASE_DAY.md). Validation log:
[`notes/phase1_validation.md`](notes/phase1_validation.md).

```powershell
.venv\Scripts\python -m pip install -e .
.venv\Scripts\python -m pytest            # fast tests
.venv\Scripts\python -m pytest -m slow    # real-data and live-archive tests
```
