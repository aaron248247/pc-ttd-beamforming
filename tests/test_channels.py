import numpy as np
import pytest

from pcttd import channels as C
from pcttd import guided as GD

F = np.round(np.arange(0.22, 0.2751, 0.0005), 6)
TAU = 3.93
F0 = 0.2495


def test_channel_lengths():
    assert C.channel_lengths(8, 1) == [8, 9, 10, 11]
    assert C.channel_lengths(4, 2, 3) == [4, 6, 8]
    with pytest.raises(ValueError):
        C.channel_lengths(8, 0)


def test_ideal_family_relative_phase_and_delay():
    H = C.ideal_channels(F, TAU)
    phi = C.relative_phase(H)
    for m in range(4):
        d = phi[m] - 2 * np.pi * F * m * TAU                 # unwrap along f leaves a constant multiple of 2 pi
        assert np.allclose(d, d[0], atol=1e-9) and abs(d[0] / (2 * np.pi) - round(d[0] / (2 * np.pi))) < 1e-9
    d = C.delays(F, phi)
    assert np.allclose(d[:, 2:-2], np.arange(4)[:, None] * TAU, atol=1e-6)
    inc = C.incremental_delays(F, H)
    assert np.allclose(inc[:, 2:-2], TAU, atol=1e-6)


def test_branch_reference_fixes_coarse_sampling():
    f = np.linspace(0.23, 0.27, 9)                         # coarse: raw unwrap along f fails (phase step > pi per sample)
    H = C.ideal_channels(f, 3.93)
    ref = 2 * np.pi * f * 3.93 * 1.01
    phi = C.relative_phase(H, ref_per_period=ref)
    assert np.allclose(phi[2], 2 * np.pi * f * 2 * 3.93, atol=1e-9)


def test_delay_error_ideal_zero_and_sign():
    H = C.ideal_channels(F, TAU)
    phi = C.relative_phase(H)
    err, rel = C.delay_error(C.delays(F, phi), TAU)
    assert np.allclose(err[:, 2:-2], 0, atol=1e-6) and np.isnan(rel[0]).all()
    H2 = C.ideal_channels(F, 1.1 * TAU)                    # 10% too much delay per step
    err2, rel2 = C.delay_error(C.delays(F, C.relative_phase(H2)), TAU)
    assert np.allclose(rel2[1:, 5:-5], 0.1, atol=1e-6) and (err2[1:, 5:-5] > 0).all()


def test_raw_phase_error_pure_intercept_and_branch():
    off = np.array([0.0, 0.3, -0.5, 2.0])
    H = C.ideal_channels(F, TAU) * np.exp(1j * off)[:, None]
    H = H / H[0] * H[0]                                    # channel 0 unchanged (offset 0)
    phi = C.relative_phase(H)
    eps = C.raw_phase_error(F, phi, TAU, F0)
    assert np.allclose(eps, off[:, None], atol=1e-9)
    # an offset of 2 pi + 0.3 is reported as 0.3 at f0
    H2 = C.ideal_channels(F, TAU) * np.exp(1j * (np.array([0, 2 * np.pi + 0.3, 0, 0])))[:, None]
    eps2 = C.raw_phase_error(F, C.relative_phase(H2) + 2 * np.pi * np.array([0, 3, 0, 0])[:, None], TAU, F0)
    assert abs(eps2[1, np.argmin(abs(F - F0))] - 0.3) < 1e-9


def test_trim_zero_at_f0_keeps_amplitude_and_delay():
    off = np.array([0.0, 0.3, -0.5, 1.2])
    amp = np.array([1.0, 0.98, 0.95, 0.9])[:, None]
    slope = 1.07                                             # delay slope differs from the target by 7%
    H = amp * C.ideal_channels(F, TAU * slope) * np.exp(1j * off)[:, None]
    phi = C.relative_phase(H)
    eps = C.raw_phase_error(F, phi, TAU, F0)
    delta = C.trim_constants(eps, F, F0)
    assert delta[0] == 0
    Ht = C.apply_trim(H, delta)
    assert np.allclose(np.abs(Ht), np.abs(H))
    phit = C.relative_phase(Ht)
    epst = C.raw_phase_error(F, phit, TAU, F0)
    i0 = np.argmin(np.abs(F - F0))
    assert np.allclose(epst[:, i0], 0, atol=1e-12)
    assert np.allclose(epst, eps - delta[:, None], atol=1e-9)
    # group delay unchanged by the trim
    assert np.allclose(C.delays(F, phit), C.delays(F, phi), atol=1e-9)
    # the residual slope error survives trimming (trimming does not hide a delay error)
    assert abs(epst[3, -1]) > 0.1


def test_trim_matches_array_trim_phase():
    from pcttd import array as A
    rng = np.random.default_rng(1)
    H = C.ideal_channels(F, TAU) * np.exp(1j * rng.uniform(-1, 1, 4))[:, None]
    eps = C.raw_phase_error(F, C.relative_phase(H), TAU, F0)
    Ht = C.apply_trim(H, C.trim_constants(eps, F, F0))
    Ha, _ = A.trim_phase(H, F, np.arange(4) * TAU, F0)
    assert np.allclose(Ht, Ha, atol=1e-12)


def test_amplitude_imbalance_is_amplitude_db():
    H = np.array([[1.0], [0.5], [1.0], [0.8]], complex)
    assert C.amplitude_imbalance_db(H)[0] == pytest.approx(20 * np.log10(2.0))   # 6.02 dB, not 3.01 dB


def test_common_mask_and_widest_run():
    a = np.array([1, 1, 1, 0, 1, 1, 1, 1, 0], bool)
    b = np.array([1, 1, 1, 1, 1, 1, 1, 0, 1], bool)
    m = C.common_usable_mask([a, b])
    f = np.arange(9) * 0.1
    assert m.tolist() == [True, True, True, False, True, True, True, False, False]
    assert C.widest_run(m, f) == (0, 2) or C.widest_run(m, f) == (4, 6)
    assert C.widest_run(np.zeros(4, bool), np.arange(4.0)) is None


def test_cascade_predict_uniform_with_interfaces():
    # lossless symmetric interface pair: S11 = r (real small), S21 = t e^{i k N}; build exact S of a length-N cavity with two identical faces
    f = np.linspace(0.23, 0.27, 7)
    k = 2 * np.pi * f * 3.9
    r1 = 0.1 * np.exp(1j * 0.3)

    def device(N):
        # two-port with face reflection r1 at both ends (lossless): exact Fabry-Perot of a face S-matrix
        t1 = np.sqrt(1 - abs(r1) ** 2)
        ph = np.exp(1j * k * N)
        # face matrix (left face): S11 = r1, S21 = t1, S22 = -conj(r1) t1/conj(t1)=-conj(r1); right face mirrored
        S11f, S21f, S22f = r1, t1, -np.conj(r1)
        d = 1 - S22f * S22f * ph * ph
        S21 = S21f * ph * S21f / d
        S11 = S11f + S21f * ph * S22f * ph * S21f / d
        return S11, S21
    a, b, c = device(8), device(9), device(11)
    S11p, S21p = C.cascade_predict(a[0], a[1], b[0], b[1], 2)
    assert np.allclose(S21p, c[1], atol=1e-9) and np.allclose(S11p, c[0], atol=1e-9)
