import numpy as np
import pandas as pd
import pytest

from companions.population import (CompletenessModel, host_completeness, mutual_inclination_posterior,
                                   occurrence_posterior)


def test_occurrence_with_full_completeness():
    p = np.ones(100)
    det = np.zeros(100, bool)
    det[:20] = True
    post = occurrence_posterior(det, p)
    assert post.median == pytest.approx(0.20, abs=0.02)
    assert 0.03 < post.sigma < 0.05


def test_occurrence_corrects_for_completeness():
    p = np.full(200, 0.5)
    det = np.zeros(200, bool)
    det[:40] = True                        # 40 of 200 at 50% completeness -> f ~ 0.4
    assert occurrence_posterior(det, p).median == pytest.approx(0.40, abs=0.05)


def test_upper_limit_without_detections():
    p = np.full(300, 0.1)                  # effective sample sum(p) = 30
    post = occurrence_posterior(np.zeros(300, bool), p)
    assert post.upper95 == pytest.approx(3.0 / 30, rel=0.15)


def test_host_completeness_limits():
    model = CompletenessModel(pd.DataFrame(dict(aligned=[True] * 4, p_lo=[0.5, 1.5, 2.5, 3.5],
                                                p_hi=[1.5, 2.5, 3.5, 4.5], n=[1] * 4,
                                                k=[3.5] * 4, ln_d50=[np.log(80)] * 4)))
    near = dict(n_obs=1000, sigma_eff=0.08, plx=50.0, mstar=1.0)      # 20 pc: everything detectable
    far = dict(n_obs=500, sigma_eff=0.3, plx=1.0, mstar=1.0)          # 1 kpc: nothing
    assert host_completeness(near, model, (1, 13), (1, 4)) > 0.95
    assert host_completeness(far, model, (1, 13), (1, 4)) < 0.01


@pytest.mark.parametrize("sigma_true", [5.0, 20.0])
def test_mutual_inclination_width_recovered(sigma_true):
    rng = np.random.default_rng(int(sigma_true))
    n = 40
    dI = np.deg2rad(rng.rayleigh(sigma_true, n))
    dO = rng.uniform(0, 2 * np.pi, n)
    i_b = np.deg2rad(rng.uniform(85, 90, n))
    cos_ic = np.cos(i_b) * np.cos(dI) + np.sin(i_b) * np.sin(dI) * np.cos(dO)
    i_c = np.rad2deg(np.arccos(cos_ic)) + rng.normal(0, 2.0, n)
    post = mutual_inclination_posterior(i_c, np.full(n, 2.0), np.rad2deg(i_b))
    lo, hi = np.percentile(post.samples, [2.5, 97.5])
    assert lo < sigma_true < hi
