"""Supercell geometry for line-defect (W1) waveguides, in lattice units (a = 1).

Pure Python/NumPy so the layout logic is unit-testable without MPB.

Conventions
-----------
* The waveguide axis is x; the supercell has period 1 along x and height ``H`` along y.
* y = 0 is the defect centre line, so the whole structure has y -> -y mirror symmetry.
* ``Scatterer(x, y, radius)``: air holes for the triangular candidate, dielectric rods for the square one.
* Positions are returned inside one supercell; periodic images are handled by the MPB solver.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

SQRT3 = math.sqrt(3.0)


@dataclass(frozen=True)
class Scatterer:
    x: float
    y: float
    radius: float


def triangular_row_pitch() -> float:
    return SQRT3 / 2.0


def triangular_w1(M: int, r: float, f1: float = 1.0, s: float = 0.0) -> tuple[float, list[Scatterer]]:
    """Triangular-lattice W1: rows j = -M..M-1 at y_j = j*sqrt(3)/2, row j offset by (j mod 2)/2 in x.

    Row j=0 is removed (the defect). Height H = M*sqrt(3) holds exactly 2M rows, which keeps the
    ABAB row stacking periodic. Rows |j|=1 get radius ``f1*r`` and are moved outward by ``s``.
    Returns (H, scatterers).
    """
    if M < 2:
        raise ValueError("M must be >= 2")
    H = M * SQRT3
    out: list[Scatterer] = []
    for j in range(-M, M):
        if j == 0:
            continue
        x = 0.5 * (j % 2)
        y = j * SQRT3 / 2.0
        rad = r
        if abs(j) == 1:
            y += math.copysign(1.0, j) * s
            rad = f1 * r
        out.append(Scatterer(x, y, rad))
    return H, out


def square_w1(M: int, r: float, f1: float = 1.0, s: float = 0.0) -> tuple[float, list[Scatterer]]:
    """Square-lattice rod W1: rows j = -M..M at y = j, row j=0 removed. Height H = 2M+1."""
    if M < 2:
        raise ValueError("M must be >= 2")
    H = 2 * M + 1.0
    out: list[Scatterer] = []
    for j in range(-M, M + 1):
        if j == 0:
            continue
        y = float(j)
        rad = r
        if abs(j) == 1:
            y += math.copysign(1.0, j) * s
            rad = f1 * r
        out.append(Scatterer(0.0, y, rad))
    return H, out


def row_pitch(lattice: str) -> float:
    return {"triangular": SQRT3 / 2.0, "square": 1.0}[lattice]


def is_mirror_symmetric(scatterers: list[Scatterer], H: float, tol: float = 1e-12) -> bool:
    """True if the structure is invariant under y -> -y modulo the supercell period H."""
    def key(sc: Scatterer):
        return (round(sc.x % 1.0, 9) % 1.0, round(sc.radius, 9))

    def wrap(y: float) -> float:
        return (y + H / 2.0) % H - H / 2.0

    for sc in scatterers:
        ym = wrap(-sc.y)
        ok = any(
            abs(wrap(o.y - ym)) < 1e-9 and abs(((o.x - sc.x) + 0.5) % 1.0 - 0.5) < tol + 1e-9
            and abs(o.radius - sc.radius) < tol
            for o in scatterers
        )
        if not ok:
            return False
    return True


def min_gap(scatterers: list[Scatterer], H: float) -> float:
    """Smallest edge-to-edge separation between any two scatterers (periodic images included).

    Negative means overlap. Used to reject unphysical sweep points.
    """
    best = math.inf
    for i, p in enumerate(scatterers):
        for j, q in enumerate(scatterers):
            if j < i:
                continue
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if i == j and dx == 0 and dy == 0:
                        continue
                    d = math.hypot(q.x + dx - p.x, q.y + dy * H - p.y)
                    best = min(best, d - p.radius - q.radius)
    return best
