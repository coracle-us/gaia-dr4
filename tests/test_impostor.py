import numpy as np
import pytest

from companions.impostor import (MainSequence, photocentre_factor, rv_signatures,
                                 twin_mimicking_signal)

AU_PER_YR_KMS = 4.740470463


@pytest.fixture(scope="module")
def ms():
    try:
        return MainSequence.load()
    except FileNotFoundError:
        pytest.skip("Pecaut & Mamajek table not downloaded")


def test_photocentre_factor_limits():
    assert photocentre_factor(1.0, 1.0) == pytest.approx(0.0)          # identical twins: no wobble
    assert photocentre_factor(0.5, 0.0) == pytest.approx(1 / 3)        # dark companion: B = q/(1+q)
    assert photocentre_factor(0.9, 0.9) == pytest.approx(0.0)          # flux ratio = mass ratio


def test_main_sequence_relation_is_sensible(ms):
    assert ms.abs_g(1.0) == pytest.approx(4.7, abs=0.3)               # Sun: M_G ~ 4.67
    m = np.linspace(0.3, 1.5, 20)
    assert np.all(np.diff(ms.abs_g(m)) < 0)                           # brighter with mass
    assert np.all(np.diff(ms.abs_rp(m)) < 0)


def test_rv_signatures_for_solar_twin_binary():
    # two 1 M_sun stars, P = 1 yr: a_rel = 2^(1/3) AU, v_rel = 2 pi a / P
    sig = rv_signatures(1.0, 1.0, 1.0, inc=np.pi / 2)
    v_rel = 2 * np.pi * 2 ** (1 / 3) * AU_PER_YR_KMS
    assert sig["k_rel"] == pytest.approx(v_rel, rel=1e-6)
    assert sig["k1"] == pytest.approx(v_rel / 2, rel=1e-6)
    assert sig["max_line_separation"] == pytest.approx(v_rel, rel=1e-6)


def test_twin_reproduces_requested_photocentre_signal(ms):
    m1, period, plx = 1.0, 3.0, 10.0
    a0_target = 0.2                                                   # mas
    tw = twin_mimicking_signal(ms, m1, period, plx, a0_target)
    assert 0.9 < tw["q"] < 1.0
    a_rel_mas = ((m1 * (1 + tw["q"])) * period**2) ** (1 / 3) * plx
    assert a_rel_mas * photocentre_factor(tw["q"], tw["flux_ratio_g"]) == pytest.approx(a0_target, rel=1e-4)
    assert tw["overluminosity_g"] > 0.6                               # near-twin: ~0.75 mag brighter
