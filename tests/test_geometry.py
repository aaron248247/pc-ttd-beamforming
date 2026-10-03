import math

import pytest

from pcttd.geometry import (SQRT3, is_mirror_symmetric, min_gap, square_w1, triangular_w1)


@pytest.mark.parametrize("M", [3, 5, 6])
def test_triangular_counts_and_height(M):
    H, sc = triangular_w1(M, 0.3)
    assert H == pytest.approx(M * SQRT3)
    assert len(sc) == 2 * M - 1          # 2M rows, one removed
    assert all(abs(s.y) > 1e-9 for s in sc)   # defect row empty


def test_triangular_symmetry_and_nn_distance():
    H, sc = triangular_w1(5, 0.3)
    assert is_mirror_symmetric(sc, H)
    # undisturbed rows: nearest-neighbour distance between adjacent rows is exactly a=1
    ys = sorted({round(s.y, 9) for s in sc})
    row = [s for s in sc if abs(s.y - SQRT3 / 2 * 2) < 1e-9][0]
    nxt = [s for s in sc if abs(s.y - SQRT3 / 2 * 3) < 1e-9][0]
    dx = abs(nxt.x - row.x)
    assert math.hypot(min(dx, 1 - dx), SQRT3 / 2) == pytest.approx(1.0)


def test_triangular_stacking_is_periodic():
    # the row that wraps (j=-M) must have the same x offset as its periodic image (j=+M)
    M = 5
    H, sc = triangular_w1(M, 0.3)
    bottom = [s for s in sc if abs(s.y + M * SQRT3 / 2) < 1e-9][0]
    assert bottom.x == pytest.approx(0.5 * (M % 2))


def test_square_counts_and_symmetry():
    H, sc = square_w1(4, 0.2)
    assert H == 9.0 and len(sc) == 8
    assert is_mirror_symmetric(sc, H)


def test_tuning_parameters_keep_symmetry_and_move_correct_rows():
    H, base = triangular_w1(5, 0.3)
    H, tuned = triangular_w1(5, 0.3, f1=1.2, s=0.1)
    assert is_mirror_symmetric(tuned, H)
    up = [s for s in tuned if s.y > 0 and s.y < 1.2][0]
    assert up.y == pytest.approx(SQRT3 / 2 + 0.1) and up.radius == pytest.approx(0.36)
    dn = [s for s in tuned if s.y < 0 and s.y > -1.2][0]
    assert dn.y == pytest.approx(-(SQRT3 / 2 + 0.1))   # outward = away from the defect on both sides
    far = [s for s in tuned if abs(s.y - SQRT3) < 1e-9][0]
    assert far.radius == pytest.approx(0.3)            # only first-neighbour rows are tuned


def test_asymmetric_structure_detected():
    H, sc = square_w1(4, 0.2)
    sc[0] = type(sc[0])(sc[0].x, sc[0].y + 0.05, sc[0].radius)
    assert not is_mirror_symmetric(sc, H)


def test_min_gap_flags_overlap():
    H, sc = triangular_w1(5, 0.3)
    assert min_gap(sc, H) == pytest.approx(1 - 0.6)
    H, sc = triangular_w1(5, 0.52)
    assert min_gap(sc, H) < 0
