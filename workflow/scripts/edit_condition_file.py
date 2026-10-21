"""Fill in the SalmonTE condition file from the 3t-seq sample sheet.

``SalmonTE.py quant`` writes ``condition.csv`` with columns ``SampleID,condition``
and every condition set to ``NA``. ``SalmonTE.py test`` is run with
``--conditions=control,treatment``, so this script sets ``condition`` to
``control`` for samples whose DESeq2 design variable equals the reference level
and ``treatment`` for all others. SalmonTE and DESeq2 thus use the same grouping.

Only the design variable column is taken from the sample sheet, so sample sheets
may carry any other columns, including one named ``condition``.
"""

import logging
from pathlib import Path

import pandas as pd

SAMPLE_ID = "SampleID"
CONDITION = "condition"
CONTROL = "control"
TREATMENT = "treatment"


def sample_sheet_index(sample_sheet: pd.DataFrame) -> pd.Series:
    """Return the sample IDs SalmonTE derives from the FASTQ files.

    Args:
        sample_sheet: The 3t-seq sample sheet.

    Returns:
        A series aligned with ``sample_sheet`` holding the SalmonTE sample ID of
        each row: ``filename`` for single-end sheets, ``filename_1`` without its
        ``_sequence`` suffix when every paired-end file has one, ``name``
        otherwise.

    Raises:
        ValueError: If the column needed to build the IDs is missing.
    """
    if "filename" in sample_sheet.columns:
        return sample_sheet["filename"]
    if "filename_1" in sample_sheet.columns and all(
        sample_sheet["filename_1"].astype(str).str.endswith("_sequence")
    ):
        return sample_sheet["filename_1"].str.split("_sequence").str[0]
    if "name" in sample_sheet.columns:
        return sample_sheet["name"]
    raise ValueError(
        "Sample sheet has neither a 'filename' nor a 'name' column; "
        f"found: {', '.join(sample_sheet.columns)}"
    )


def edit_condition(
    condition_sheet: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    variable: str,
    reference_level: str,
) -> pd.DataFrame:
    """Build the SalmonTE condition table.

    Args:
        condition_sheet: ``condition.csv`` as written by ``SalmonTE.py quant``.
        sample_sheet: The 3t-seq sample sheet.
        variable: Sample sheet column holding the DESeq2 design variable.
        reference_level: Level of ``variable`` mapped to ``control``.

    Returns:
        A ``SampleID,condition`` table in the row order of ``condition_sheet``.

    Raises:
        ValueError: If ``variable`` is not a sample sheet column, sample IDs
            are duplicated, a SalmonTE sample is missing from the sample sheet,
            or no sample has ``reference_level``.
    """
    if SAMPLE_ID not in condition_sheet.columns:
        raise ValueError(
            f"SalmonTE condition file has no '{SAMPLE_ID}' column; "
            f"found: {', '.join(condition_sheet.columns)}"
        )
    if variable not in sample_sheet.columns:
        raise ValueError(
            f"DESeq2 design variable '{variable}' is not a column of the sample "
            f"sheet (columns: {', '.join(sample_sheet.columns)}). Set "
            "deseq2.variable to the sample sheet column to compare."
        )

    levels = pd.Series(
        sample_sheet[variable].to_numpy(),
        index=sample_sheet_index(sample_sheet).to_numpy(),
    )
    duplicated = levels.index[levels.index.duplicated()].unique()
    if not duplicated.empty:
        raise ValueError(
            f"Duplicated sample IDs in the sample sheet: {', '.join(duplicated)}"
        )
    sample_ids = condition_sheet[SAMPLE_ID]
    missing = sample_ids[~sample_ids.isin(levels.index)]
    if not missing.empty:
        raise ValueError(
            f"SalmonTE samples not found in the sample sheet: {', '.join(missing)}"
        )

    levels = levels.loc[sample_ids].astype(str)
    is_control = levels == str(reference_level)
    if not is_control.any():
        raise ValueError(
            f"No sample has reference level '{reference_level}' in column "
            f"'{variable}' (levels: {', '.join(sorted(levels.unique()))})"
        )

    return pd.DataFrame(
        {
            SAMPLE_ID: sample_ids.to_numpy(),
            CONDITION: [CONTROL if c else TREATMENT for c in is_control],
        }
    )


def main(
    condition_file: Path, sample_sheet_file: Path, variable: str, reference_level: str
) -> None:
    """Rewrite ``condition_file`` in place with control/treatment labels."""
    condition_sheet = pd.read_csv(condition_file)
    sample_sheet = pd.read_csv(sample_sheet_file)
    out = edit_condition(condition_sheet, sample_sheet, variable, reference_level)
    logging.info("Writing %s:\n%s", condition_file, out.to_string(index=False))
    out.to_csv(condition_file, index=False)


if "snakemake" in globals():
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    try:
        main(
            Path(snakemake.input.condition_file),
            Path(snakemake.input.sample_sheet),
            snakemake.params.variable,
            snakemake.params.reference_level,
        )
    except ValueError:
        logging.exception("Could not write the SalmonTE condition file")
        raise
