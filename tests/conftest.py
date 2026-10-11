from pathlib import Path

import pytest
from astropy.io.votable import parse_single_table

ROOT = Path(__file__).resolve().parents[1]
PRERELEASE_XML = (ROOT / "data/raw/dr4_prerelease/epoch_astrometry/"
                  "GAIA_DR4_PRERELEASE_EPOCH_ASTROMETRY_RAW.xml")
GAIA4 = 1457486023639239296
BH3 = 4318465066420528000
HD183633 = 4181040337841125632


@pytest.fixture(scope="session")
def prerelease_table():
    if not PRERELEASE_XML.exists():
        pytest.skip("DR4 pre-release sample not downloaded (see notes/dr4_data_model.md)")
    return parse_single_table(PRERELEASE_XML).to_table()


@pytest.fixture(scope="session")
def gaiasupdate_fit(prerelease_table):
    """ESA gaiasupdate DR4-like (6p_constrained_colour) fit, keyed by source_id."""
    from gaiasupdate.epoch_astrometry import GaiaEpochAstrometryArchive
    cache = {}

    def fit(source_id):
        if source_id not in cache:
            df = prerelease_table[prerelease_table["source_id"] == source_id].to_pandas()
            cache[source_id] = GaiaEpochAstrometryArchive.supdate(df, source_id)
        return cache[source_id]
    return fit
