import pandas as pd
from pathlib import Path

serie_folder = Path(snakemake.input[0])
sample_sheet_path = Path(snakemake.input[1])

experimental_variable = snakemake.params.variable
reference_level = snakemake.params.reference_level

salmon_condition_sheet = pd.read_csv(serie_folder.joinpath("condition.csv"))
salmon_condition_sheet = salmon_condition_sheet.set_index("SampleID")

sample_sheet = pd.read_csv(sample_sheet_path)

if "filename" in sample_sheet.columns:
    sample_sheet = sample_sheet.set_index("filename")
else:
    if all(sample_sheet.filename_1.apply(lambda x: x.endswith("_sequence"))):
        sample_sheet["tmp"] = sample_sheet.filename_1.apply(
            lambda x: x.split("_sequence")[0]
        )
        sample_sheet = sample_sheet.set_index("tmp")
    else:
        sample_sheet = sample_sheet.set_index("name")

sample_sheet = sample_sheet.reindex(index=salmon_condition_sheet.index)

joined = sample_sheet.join(salmon_condition_sheet)

joined["condition"] = joined.apply(
    lambda row: "control"
    if row[snakemake.params.variable] == reference_level
    else "treatment",
    axis=1,
)
joined = joined.reset_index()

out = joined[["SampleID", "condition"]]
out.to_csv(serie_folder.joinpath("condition.csv"), index=False)
