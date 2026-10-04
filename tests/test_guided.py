"""Pure-NumPy tests of the guided-mode projection, 2-port algebra and delay estimators (``pcttd.guided``)."""
import numpy as np
import pytest

from pcttd import extract as X
from pcttd import guided as GD

RNG = np.random.default_rng(7)


def _mode(ny=64):
    y = (np.arange(ny) / ny - 0.5) * 6.0
    mEy = (np.exp(-y**2) * (1 + 0.3j)) + 0.1j * y * np.exp(-y**2 / 2)
    mHz = np.exp(-(y - 0.2)**2 / 1.5) * (0.8 - 0.2j)
    return y, mEy, mHz


def test_projection_recovers_known_forward_and_backward_amplitudes_exactly():
    y, mEy, mHz = _mode()
    w = np.full(len(y), y[1] - y[0])
    af0, ab0 = 0.7 - 0.4j, -0.2 + 0.9j
    Ey = af0 * mEy + ab0 * np.conj(mEy)
    Hz = af0 * mHz + ab0 * (-np.conj(mHz))
    af, ab, res = GD.project(Ey, Hz, mEy, mHz, w)
    assert np.allclose(af, af0) and np.allclose(ab, ab0)
    assert res[0] < 1e-12


def test_projection_residual_detects_a_foreign_field_and_is_scale_invariant():
    y, mEy, mHz = _mode()
    w = np.full(len(y), y[1] - y[0])
    other = np.cos(3 * y) * np.exp(-y**2 / 4)
    Ey = mEy + 0.3 * other
    Hz = mHz + 0.3 * other
    _, _, r1 = GD.project(Ey, Hz, mEy, mHz, w)
    _, _, r2 = GD.project(5 * Ey, 5 * Hz, mEy, mHz, w)
    assert 1e-3 < r1[0] < 1.0
    assert np.allclose(r1, r2)


def test_forward_and_backward_modes_are_orthogonal_with_opposite_norm():
    y, mEy, mHz = _mode()
    w = np.full(len(y), 1.0)
    # N_f from projecting the forward mode onto itself, N_b from the partner
    bEy, bHz = np.conj(mEy), -np.conj(mHz)
    Nf = np.sum(w * (np.conj(mEy) * mHz + mEy * np.conj(mHz)))
    Nb = np.sum(w * (np.conj(bEy) * bHz + bEy * np.conj(bHz)))
    assert np.isclose(Nb, -Nf)
    Of_b = np.sum(w * (np.conj(mEy) * bHz + bEy * np.conj(mHz)))
    assert abs(Of_b) < 1e-12


def test_resample_periodic_reproduces_a_periodic_function():
    H, n = 5.0, 80
    y = (np.arange(n) / n - 0.5) * H
    f = lambda t: np.sin(2 * np.pi * t / H) + 0.5j * np.cos(4 * np.pi * t / H)
    yn = np.linspace(-3, 3, 37)                     # includes points outside [-H/2, H/2)
    assert np.allclose(GD.resample_periodic(y, f(y), yn, H), f(yn), atol=2e-4)


def _lossless_2port(rng):
    """Random lossless reciprocal-symmetric-phase 2-port: S11, S21 with S22 from unitarity."""
    R = rng.uniform(0.05, 0.6)
    S11 = np.sqrt(R) * np.exp(1j * rng.uniform(-np.pi, np.pi))
    S21 = np.sqrt(1 - R) * np.exp(1j * rng.uniform(-np.pi, np.pi))
    return S11, S21, GD.s22_from_unitarity(S11, S21)


def test_s22_from_unitarity_gives_a_unitary_matrix():
    S11, S21, S22 = _lossless_2port(RNG)
    S = np.array([[S11, S21], [S21, S22]])
    assert np.allclose(S.conj().T @ S, np.eye(2))


def test_t_matrix_roundtrip_and_cascade_rule():
    A = _lossless_2port(RNG); B = _lossless_2port(RNG)
    TA, TB = GD.t_matrix(*A), GD.t_matrix(*B)
    s11, s21, s12, s22 = GD.s_from_t(TA)
    assert np.allclose([s11, s21, s12, s22], [A[0], A[1], A[1], A[2]])
    # cascade A then B, by the explicit Redheffer-type formulas
    d = 1 - A[2] * B[0]
    S21c = A[1] * B[1] / d
    S11c = A[0] + A[1] * A[1] * B[0] / d
    s11c, s21c, _, _ = GD.s_from_t(TB @ TA)
    assert np.isclose(s21c, S21c) and np.isclose(s11c, S11c)


def test_propagate_t_is_a_pure_phase_two_port_with_meep_sign():
    k, L = 0.37, 5.0
    s11, s21, s12, s22 = GD.s_from_t(GD.propagate_T(k, L))
    assert np.isclose(s21, np.exp(1j * k * L)) and abs(s11) < 1e-15 and abs(s22) < 1e-15


def test_section_phase_from_lengths_is_independent_of_the_interfaces():
    S11, S21, S22 = _lossless_2port(RNG)
    Tr = GD.t_matrix(S11, S21, S22)                       # interface (device-lead step), same on both ends
    Tl = np.linalg.inv(Tr)
    k = np.array([0.31, 1.7, -0.4])
    def device(N):
        return np.array([Tr @ GD.propagate_T(kk, N) @ Tl for kk in k])
    th, mod = GD.section_phase_from_lengths(device(9), device(5))
    want = np.angle(np.exp(1j * k * 4))
    assert np.allclose(th, want, atol=1e-9) and np.allclose(mod, 1.0)


def test_local_poly_derivative_is_exact_for_quadratics_including_ends():
    x = np.linspace(0, 2, 41)
    y = 3 + 2 * x - 1.5 * x**2
    assert np.allclose(GD.local_poly_derivative(x, y), 2 - 3 * x, atol=1e-9)


@pytest.mark.parametrize("method", ["central", "poly"])
def test_delay_from_phase_recovers_positive_tau_for_meep_convention(method):
    f = np.linspace(0.2, 0.3, 201)
    tau = 37.5
    ph = X.omega(f) * tau + 0.9                   # arg H of exp(+i omega tau)
    t = GD.delay_from_phase(f, ph, method=method)
    assert np.allclose(t[4:-4], tau, atol=1e-8)


def test_unwrap_theta_with_reference_fixes_branch_on_coarse_grids():
    true = np.linspace(0, 40, 9)                  # steps of ~5 rad: plain unwrap would fail
    wrapped = np.angle(np.exp(1j * true))
    assert not np.allclose(np.unwrap(wrapped), true)
    assert np.allclose(GD.unwrap_theta(wrapped, ref=true + 0.4), true)
    assert np.allclose(GD.unwrap_theta(np.angle(np.exp(1j * np.linspace(0, 20, 400)))), np.linspace(0, 20, 400))
