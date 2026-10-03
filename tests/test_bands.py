import numpy as np
import pytest

from pcttd import bands as B


def make_band(w_of_k, k):
    w = w_of_k(k)
    h = 1e-6
    vg = (w_of_k(k + h) - w_of_k(k - h)) / (2 * h)
    return w, vg


def build(bands, k):
    """bands: list of (w_of_k, parity). Returns branch dicts as class_branches would."""
    br = []
    for f, p in bands:
        w, vg = make_band(f, k)
        br.append(dict(parity=p, w=w, vg=vg, strip=np.ones_like(k) * 0.9, edge=np.zeros_like(k)))
    return br


K = np.linspace(0.0, 0.5, 41)


def test_find_gaps():
    f = np.array([[0.1, 0.5, 0.7], [0.2, 0.6, 0.8], [0.15, 0.55, 0.75]])
    g = B.find_gaps(f)
    assert len(g) == 2
    assert g[0].n == 1 and g[0].w_lo == pytest.approx(0.2) and g[0].w_hi == pytest.approx(0.5)
    assert g[0].rel == pytest.approx(0.3 / 0.35)
    # overlapping bands -> no gap
    assert B.find_gaps(np.array([[0.1, 0.3], [0.4, 0.5]])) == []


def test_hermite_reproduces_smooth_function_and_slope():
    f = lambda k: 0.25 + 0.1 * np.sin(np.pi * k)
    w, vg = make_band(f, K)
    kf, wf, df = B.hermite_resample(K, w, vg, 24)
    assert np.max(np.abs(wf - f(kf))) < 1e-8
    assert np.max(np.abs(df - 0.1 * np.pi * np.cos(np.pi * kf))) < 1e-5


def test_constant_group_index_gbp_matches_analytic():
    n_g = 5.0
    f = lambda k: 0.30 + k / n_g          # linear band, n_g = 5, spans 0.30..0.40
    br = build([(f, +1)], K)
    levels = np.linspace(0.301, 0.399, 400)
    c = B.find_crossings(K, br, levels)
    assert (c.counts() == 1).all()
    assert np.allclose(1 / c.vg, n_g, rtol=1e-4)
    strict_any, strict_even, even_only = B.single_mode_masks(c)
    assert strict_even.all()
    vg, _, _ = B.per_level_single(c, strict_even[c.level_idx])
    res = B.best_flat_window(levels, 1 / vg, n_min=5 - 1e-3)
    lo, hi = levels[0], levels[-1]
    assert res["gbp"] == pytest.approx(n_g * (hi - lo) / (0.5 * (hi + lo)), rel=1e-3)
    assert res["n_mean"] == pytest.approx(n_g, rel=1e-4)
    assert B.best_flat_window(levels, 1 / vg, n_min=6.0) is None      # below the n_min threshold


def test_flat_window_rejects_non_flat_band():
    # n_g varies 2 -> 20 across the window: only short sub-windows can satisfy +-10%
    w = np.linspace(0.3, 0.4, 500)
    ng = np.linspace(2, 20, 500)
    r = B.best_flat_window(w, ng, n_min=0.0)
    # analytic bound: a +-10% window around n_mean <= 20 spans at most dn = 0.2*20 = 4 of the 18-wide ramp,
    # i.e. dw <= (4/18)*0.1 -> frac_bw <= (4/18)*0.1/0.35 ~ 0.0635 (full window would be 0.286)
    assert r is not None and 0 < r["frac_bw"] <= (4 / 18) * 0.1 / 0.35 + 1e-9
    # every point inside the returned window obeys the tolerance
    i = np.argmin(np.abs(w - r["w_lo"])); j = np.argmin(np.abs(w - r["w_hi"]))
    seg = ng[i:j + 1]
    assert seg.max() <= 1.1 * seg.mean() + 1e-9 and seg.min() >= 0.9 * seg.mean() - 1e-9


def test_single_mode_windows_and_parity_split():
    A = lambda k: 0.25 + 0.10 * k / 0.5            # even, covers 0.25..0.35
    Bo = lambda k: 0.30 + 0.03 * k / 0.5           # odd, covers 0.30..0.33
    br = build([(A, +1), (Bo, -1)], K)
    levels = np.linspace(0.24, 0.36, 1201)
    c = B.find_crossings(K, br, levels)
    strict_any, strict_even, even_only = B.single_mode_masks(c)
    n = c.counts()
    # two modes where both exist
    assert (n[(levels > 0.301) & (levels < 0.329)] == 2).all()
    assert (n[(levels > 0.251) & (levels < 0.299)] == 1).all()
    assert not strict_any[(levels > 0.301) & (levels < 0.329)].any()
    assert strict_even[(levels > 0.251) & (levels < 0.299)].all()
    # odd band does not break the even-only window
    assert even_only[(levels > 0.251) & (levels < 0.349)].all()
    r = B.runs(strict_even)
    assert len(r) == 2                      # [0.25,0.30) and (0.33,0.35]
    w1 = B.frac_width(levels, *r[0])
    assert w1 == pytest.approx(0.05 / 0.275, rel=0.05)


def test_negative_slope_in_folded_zone_is_a_valid_forward_mode():
    # Revision 1: dw/dk < 0 on [0, 0.5] is a positive-power-flow mode with true wavevector in (0.5, 1).
    # Its time-reversed partner exists, so one crossing = one forward mode; n_g = 1/|v_g|.
    back = lambda k: 0.35 - 0.1 * k / 0.5
    br = build([(back, +1)], K)
    levels = np.linspace(0.26, 0.34, 100)
    c = B.find_crossings(K, br, levels)
    strict_any, strict_even, even_only = B.single_mode_masks(c)
    assert (c.counts() == 1).all()
    assert strict_any.all() and strict_even.all() and even_only.all()
    assert np.allclose(1 / np.abs(c.vg), 5.0, rtol=1e-4)


def test_turning_point_gives_two_modes_not_one():
    # band with a maximum inside the zone: levels just below the maximum are crossed twice
    f = lambda k: 0.30 - 0.4 * (k - 0.25) ** 2
    br = build([(f, +1)], K)
    levels = np.array([0.2999, 0.28])
    c = B.find_crossings(K, br, levels)
    assert list(c.counts()) == [2, 2]
    sa, se, eo = B.single_mode_masks(c)
    assert not sa.any()


def test_zone_edge_slow_light_group_index_diverges():
    f = lambda k: 0.25 + 0.1 * np.sin(np.pi * k)
    br = build([(f, +1)], K)
    levels = np.array([0.30, 0.34, 0.349])
    c = B.find_crossings(K, br, levels)
    ng = 1 / c.vg
    expected = 1.0 / (0.1 * np.pi * np.cos(np.arcsin((levels - 0.25) / 0.1)))
    assert np.allclose(ng, expected, rtol=1e-3)
    assert ng[2] > ng[1] > ng[0]


def test_class_branches_separates_parity_and_avoids_kink():
    # even and odd bands that cross exactly: per-parity sorting keeps each straight
    nk = len(K)
    fa = 0.30 + 0.1 * K            # even, rising
    fb = 0.40 - 0.1 * K            # odd, falling; crosses fa at K=0.5*... fa=fb at K=0.5
    freq = np.stack([fa, fb, np.full(nk, 0.9)], axis=1)
    vg = np.stack([np.full(nk, 0.1), np.full(nk, -0.1), np.zeros(nk)], axis=1)
    par = np.tile([+1, -1, +1], (nk, 1))
    z = np.zeros_like(freq)
    br = B.class_branches(freq, vg, par, z, z)
    even = [b for b in br if b["parity"] == +1]
    assert len(even) >= 1
    assert np.allclose(even[0]["w"], fa) and np.allclose(even[0]["vg"], 0.1)
