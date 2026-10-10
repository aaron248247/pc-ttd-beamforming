"""Integrity of the Phase 4 -> Phase 5 handoff dataset (skipped if the Phase 4 results are absent)."""
import json
from pathlib import Path

import numpy as np
import pytest

from pcttd import channels as C

P = Path(__file__).resolve().parents[1] / "results" / "phase4" / "phase4_handoff.npz"
pytestmark = pytest.mark.skipif(not P.exists(), reason="Phase 4 handoff not generated")


@pytest.fixture(scope="module")
def d():
    z = np.load(P)
    return {k: z[k] for k in z.files}


def test_shapes_and_contents(d):
    f = d["f"]
    assert d["H_raw"].shape == (4, len(f)) == d["H_trim"].shape
    assert list(d["N"]) == [8, 9, 10, 11]
    for k in ("amp_raw", "phase_raw_wrapped", "phase_abs_unwrapped", "phase_rel_unwrapped", "mask_common", "mask_validated_domain", "f0", "dtau_target", "S11_raw", "provenance_json"):
        assert k in d
    prov = json.loads(str(d["provenance_json"]))
    assert "exp(+i" in prov["sign_convention"] and prov["N"] == [8, 9, 10, 11]


def test_amplitude_and_phase_views_match_H(d):
    assert np.allclose(d["amp_raw"], np.abs(d["H_raw"]))
    assert np.allclose(np.exp(1j * d["phase_raw_wrapped"]), d["H_raw"] / np.abs(d["H_raw"]))
    assert np.allclose(np.exp(1j * d["phase_abs_unwrapped"]), d["H_raw"] / np.abs(d["H_raw"]))


def test_trim_is_one_constant_per_channel(d):
    assert np.allclose(d["H_trim"], d["H_raw"] * np.exp(-1j * d["trim_delta"])[:, None])
    assert np.allclose(np.abs(d["H_trim"]), np.abs(d["H_raw"]))
    assert d["trim_delta"][0] == 0


def test_trim_zero_error_at_f0_and_sign(d):
    f = d["f"]
    phi = C.relative_phase(d["H_trim"])
    eps = C.raw_phase_error(f, phi, float(d["dtau_target"]), float(d["f0"]))
    i0 = int(np.argmin(abs(f - float(d["f0"]))))
    assert np.allclose(eps[:, i0], 0, atol=1e-9)
    sel = d["mask_common"].astype(bool)
    assert np.abs(eps[:, sel]).max() < 0.2                       # trimmed channels track the ideal TTD closely
    assert (np.gradient(d["phase_rel_unwrapped"][1], 2 * np.pi * f)[sel] > 0).all()   # tau_g = +d arg H / d omega > 0 (MEEP_SIGN = +1)


def test_masks_nested(d):
    assert not (d["mask_common"] & ~d["mask_validated_domain"]).any()
    assert d["mask_common"].any()
