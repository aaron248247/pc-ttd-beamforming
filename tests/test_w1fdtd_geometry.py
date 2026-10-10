"""The FDTD device-section hole layout must reproduce the MPB supercell geometry (``geometry.triangular_w1``) for the engineered rows.
Regression test for a sign bug (copysign(s, j) discarded the sign of s) that made negative s an outward shift in the Phase 3 first-pass runs."""
import pytest

from pcttd import geometry as G
from pcttd import w1fdtd as W


@pytest.mark.parametrize("f1,s", [(1.0, 0.0), (1.2, -0.05), (1.2, -0.10), (1.1, 0.05)])
def test_device_rows_match_mpb_geometry(f1, s):
    cfg = W.W1Cfg(N=4, f1=f1, s=s)
    _, sc = G.triangular_w1(6, cfg.r, f1, s)
    mpb_rows = {(j_sign, round(x % 1.0, 9)): (y, rad) for (x, y, rad) in [(a.x, a.y, a.radius) for a in sc] for j_sign in [1 if y > 0 else -1] if abs(y) < 1.2}
    dev = [(x, y, rad) for (x, y, rad) in W.hole_list(cfg) if abs(y) < 1.2 and 0.0 <= x < cfg.N]
    assert len(dev) == 2 * cfg.N
    for x, y, rad in dev:
        y_ref, rad_ref = mpb_rows[(1 if y > 0 else -1, round(x % 1.0, 9))]
        assert y == pytest.approx(y_ref) and rad == pytest.approx(rad_ref)
    lead = [(x, y, rad) for (x, y, rad) in W.hole_list(cfg) if abs(y) < 1.2 and x < -0.1]
    assert all(rad == pytest.approx(cfg.r) and abs(abs(y) - 3 ** 0.5 / 2) < 1e-12 for x, y, rad in lead)
