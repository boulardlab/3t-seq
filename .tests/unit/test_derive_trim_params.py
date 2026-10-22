"""Tests for workflow/scripts/derive_trim_params.py FastQC parsing (issue #51)."""

import importlib.util
import sys
import zipfile
from pathlib import Path
from types import ModuleType

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow/scripts"

FASTQC_DATA = """##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Sequence length\t35-150
>>END_MODULE
"""


def load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dtp = load("derive_trim_params")


def test_zip_renamed_by_the_workflow(tmp_path: Path) -> None:
    # The wrapper renames the zip; the folder inside keeps the FASTQ name.
    zip_path = tmp_path / "S1_raw_1_fastqc.zip"
    with zipfile.ZipFile(zip_path, "w") as z:
        z.writestr("S1_1_fastqc/fastqc_data.txt", FASTQC_DATA)
    assert dtp.parse_fastqc_data(str(zip_path))["max_length"] == 150


def test_zip_without_fastqc_data_fails(tmp_path: Path) -> None:
    zip_path = tmp_path / "empty_fastqc.zip"
    with zipfile.ZipFile(zip_path, "w") as z:
        z.writestr("x/summary.txt", "")
    with pytest.raises(ValueError):
        dtp.parse_fastqc_data(str(zip_path))
