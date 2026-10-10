"""Phase 5: consistency of the Phase 4 handoff as used by the array analysis, and of the saved Phase 5 results (skipped if absent)."""
import json
from pathlib import Path

import numpy as np
import pytest

R = Path(__file__).resolve().parents[1] / "results"
P4 = R / "phase4" / "phase4_handoff.npz"
P5 = R / "phase5" / "5a_array.npz"


@pytest.fixture(scope="module")
def d():
    z = np.load(P4)
    return {k: z[k] for k in z.files}


@pytest.mark.skipif(not P4.exists(), reason="Phase 4 handoff absent")
class TestHandoff:
    def test_dimensions_and_grid(self, d):
        f = d["f"]
        assert d["H_raw"].shape == d["H_trim"].shape == d["H_raw_res40"].shape == (4, len(f))
        assert np.all(np.diff(f) > 0) and np.allclose(np.diff(f), 5e-4)
        assert d["mask_common"].sum() > 50 and f[d["mask_common"]].min() >= 0.231 - 1e-9 and f[d["mask_common"]].max() <= 0.268 + 1e-9
        assert float(d["f0"]) == 0.2495 and float(d["dtau_target"]) == 3.9323

    @pytest.mark.parametrize("sfx", ["", "_res40"])
    def test_trim_is_constant_phase_per_channel(self, d, sfx):
        ratio = d["H_trim" + sfx] / d["H_raw" + sfx]
        assert np.allclose(np.abs(ratio), 1.0, atol=1e-12)
        assert np.allclose(ratio, ratio[:, [0]], atol=1e-12)                   # frequency-independent

    @pytest.mark.parametrize("sfx", ["", "_res40"])
    def test_trim_preserves_amplitude_and_group_delay(self, d, sfx):
        Hr, Ht, f = d["H_raw" + sfx], d["H_trim" + sfx], d["f"]
        assert np.allclose(np.abs(Hr), np.abs(Ht), atol=1e-12)
        slope = lambda H: np.gradient(np.unwrap(np.angle(H), axis=1), f, axis=1)
        assert np.allclose(slope(Hr), slope(Ht), atol=1e-8)

    def test_trim_constant_matches_handoff_delta_res30(self, d):
        ratio = (d["H_trim"] / d["H_raw"])[:, 0]
        assert np.allclose(ratio, np.exp(-1j * d["trim_delta"]), atol=1e-12)


@pytest.mark.skipif(not P5.exists(), reason="Phase 5 results not generated")
def test_phase5_results_reload_and_analytic_checks():
    z = np.load(P5)
    j = json.loads((P5.with_suffix(".json")).read_text())
    assert j["fixed_analytic_check_max_dev_deg"] < 0.01
    assert np.abs(z["res30_ideal_theta_peak"] - 30.0).max() < 1e-3
    assert z["patterns_res30"].shape[:2] == (4, 3)
    assert j["setup"]["d"] == pytest.approx(7.8646)
