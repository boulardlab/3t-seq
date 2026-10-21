"""Tests for workflow/scripts/trimmomatic_steps.py (issue #48)."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow/scripts"


def load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ts = load("trimmomatic_steps")
dtp = load("derive_trim_params")

CLIP_PE = "ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10"


@pytest.mark.parametrize("paired,suffix", [(True, "PE"), (False, "SE")])
def test_default_without_extra_params_is_unchanged(paired, suffix):
    # Same string as the pre-#48 default fallback in common.smk.
    assert ts.merge_steps(ts.default_steps(paired), "") == (
        f"ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-{suffix}.fa"
        ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    )


def test_user_step_replaces_base_step_of_same_type():
    out = ts.merge_steps(ts.default_steps(True), "SLIDINGWINDOW:4:25").split()
    assert "SLIDINGWINDOW:4:25" in out
    assert "SLIDINGWINDOW:4:15" not in out


def test_order_clip_first_minlen_last_user_steps_after_base():
    out = ts.merge_steps(ts.default_steps(True), "CROP:50 HEADCROP:5").split()
    assert out[0] == CLIP_PE
    assert out[-1].startswith("MINLEN:")
    assert out.index("SLIDINGWINDOW:4:15") < out.index("CROP:50")
    assert out.index("CROP:50") < out.index("HEADCROP:5")


def test_user_illuminaclip_and_minlen_keep_their_slots():
    out = ts.merge_steps(
        ts.default_steps(True), "MINLEN:20 CROP:50 ILLUMINACLIP:custom.fa:2:30:10"
    ).split()
    assert out[0] == "ILLUMINACLIP:custom.fa:2:30:10"
    assert out[-1] == "MINLEN:20"


def test_base_minlen_is_capped_by_crop():
    base = ["ILLUMINACLIP:x.fa:2:30:10", "SLIDINGWINDOW:4:20", "MINLEN:52"]
    assert ts.merge_steps(base, "CROP:50").split()[-1] == "MINLEN:50"


def test_user_minlen_is_not_capped():
    base = ["ILLUMINACLIP:x.fa:2:30:10", "MINLEN:36"]
    assert ts.merge_steps(base, "CROP:50 MINLEN:60").split()[-1] == "MINLEN:60"


def test_adaptive_150bp_with_crop_keeps_reads():
    # 150 bp reads: derived MINLEN would be 52 > CROP 50 and drop every read.
    fastqc = [{"max_length": 150, "quality_profile": [{"p10": 30}]}]
    out = dtp.derive_params(fastqc, is_paired=True, extra_params="CROP:50").split()
    assert out[-1] == "MINLEN:50"
    assert out[0].startswith("ILLUMINACLIP:")
    assert "CROP:50" in out
