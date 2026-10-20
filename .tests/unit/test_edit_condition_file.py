"""Regression tests for workflow/scripts/edit_condition_file.py (issue #8)."""

import importlib.util
from io import StringIO
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "workflow/scripts/edit_condition_file.py"
)


def load_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("edit_condition_file", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ecf = load_script()

# condition.csv as written by `SalmonTE.py quant`.
SALMONTE_CONDITION = pd.DataFrame(
    {"SampleID": ["KO_1", "WT_1", "KO_2", "WT_2"], "condition": ["NA"] * 4}
)
EXPECTED = pd.DataFrame(
    {
        "SampleID": ["KO_1", "WT_1", "KO_2", "WT_2"],
        "condition": ["treatment", "control", "treatment", "control"],
    }
)


def sheet(csv: str) -> pd.DataFrame:
    return pd.read_csv(StringIO(csv))


def test_without_condition_column() -> None:
    sample_sheet = sheet(
        "name,filename_1,filename_2,genotype\n"
        "WT_1,a_1,a_2,WT\nWT_2,b_1,b_2,WT\nKO_1,c_1,c_2,KO\nKO_2,d_1,d_2,KO\n"
    )
    out = ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")
    pd.testing.assert_frame_equal(out, EXPECTED)


def test_condition_column_is_design_variable() -> None:
    # The sample sheet from issue #8.
    sample_sheet = sheet(
        "name,filename_1,filename_2,condition\n"
        "WT_1,a_1,a_2,CTRL\nKO_1,c_1,c_2,KD\nWT_2,b_1,b_2,CTRL\nKO_2,d_1,d_2,KD\n"
    )
    out = ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "condition", "CTRL")
    pd.testing.assert_frame_equal(out, EXPECTED)


def test_condition_column_is_not_design_variable() -> None:
    # `condition` must not leak into SalmonTE's grouping: deseq2.variable wins.
    sample_sheet = sheet(
        "name,filename_1,filename_2,genotype,condition\n"
        "WT_1,a_1,a_2,WT,KO\nWT_2,b_1,b_2,WT,KO\nKO_1,c_1,c_2,KO,WT\nKO_2,d_1,d_2,KO,WT\n"
    )
    out = ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")
    pd.testing.assert_frame_equal(out, EXPECTED)


def test_extra_shared_columns() -> None:
    # Every SalmonTE column name also appears in the sample sheet.
    sample_sheet = sheet(
        "SampleID,filename,condition,genotype\n"
        "x1,KO_1,a,KO\nx2,WT_1,b,WT\nx3,KO_2,c,KO\nx4,WT_2,d,WT\n"
    )
    out = ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")
    pd.testing.assert_frame_equal(out, EXPECTED)


def test_sequence_suffix_ids() -> None:
    sample_sheet = sheet(
        "name,filename_1,filename_2,genotype\n"
        "s1,WT_1_sequence,WT_1_r2,WT\ns2,WT_2_sequence,WT_2_r2,WT\n"
        "s3,KO_1_sequence,KO_1_r2,KO\ns4,KO_2_sequence,KO_2_r2,KO\n"
    )
    out = ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")
    pd.testing.assert_frame_equal(out, EXPECTED)


def test_missing_design_variable() -> None:
    sample_sheet = sheet("name,filename_1,filename_2,condition\nWT_1,a,b,CTRL\n")
    with pytest.raises(ValueError, match="design variable 'genotype'"):
        ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")


def test_sample_missing_from_sample_sheet() -> None:
    sample_sheet = sheet(
        "name,filename_1,filename_2,genotype\n"
        "WT_1,a_1,a_2,WT\nWT_2,b_1,b_2,WT\nKO_1,c_1,c_2,KO\n"
    )
    with pytest.raises(ValueError, match="not found in the sample sheet: KO_2"):
        ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "WT")


def test_reference_level_absent() -> None:
    sample_sheet = sheet(
        "name,filename_1,filename_2,genotype\n"
        "WT_1,a_1,a_2,WT\nWT_2,b_1,b_2,WT\nKO_1,c_1,c_2,KO\nKO_2,d_1,d_2,KO\n"
    )
    with pytest.raises(ValueError, match="reference level 'ctrl'"):
        ecf.edit_condition(SALMONTE_CONDITION, sample_sheet, "genotype", "ctrl")


def test_main_rewrites_file(tmp_path: Path) -> None:
    condition_file = tmp_path / "condition.csv"
    SALMONTE_CONDITION.to_csv(condition_file, index=False)
    sample_sheet_file = tmp_path / "samples.csv"
    sample_sheet_file.write_text(
        "name,filename_1,filename_2,condition\n"
        "WT_1,a_1,a_2,CTRL\nWT_2,b_1,b_2,CTRL\nKO_1,c_1,c_2,KD\nKO_2,d_1,d_2,KD\n"
    )
    ecf.main(condition_file, sample_sheet_file, "condition", "CTRL")
    pd.testing.assert_frame_equal(pd.read_csv(condition_file), EXPECTED)
