import os
import json
import itertools
from yte import process_yaml


def get_name(x: str) -> str:
    bn = os.path.basename(x)
    name = os.path.splitext(bn)[0]
    return name


datasets = []
views = []
max_rows = 0
for dataset_path in snakemake.input["datasets"]:
    dataset_name = get_name(str(dataset_path))
    dataset = {"name": dataset_name, "path": dataset_path}
    datasets.append(dataset)

    view = {"name": dataset_name, "dataset": dataset_name, "type": "table"}
    views.append(view)

    with open(dataset_path) as fh:
        nlines = len(fh.readlines())
        if nlines > max_rows:
            max_rows = nlines

for dataset_path, view_spec in itertools.product(
    snakemake.input["datasets"], snakemake.params["view_specs"]
):
    dataset_name = get_name(str(dataset_path))
    spec_name = get_name(str(view_spec))
    view_name = "-".join([dataset_name, spec_name])
    view = {
        "name": view_name,
        "dataset": dataset_name,
        "spec_path": view_spec,
        "type": "plot",
    }
    views.append(view)

variables = {
    "plot_name": snakemake.params["plot_name"],
    "default_view": views[0]["name"],
    "datasets": datasets,
    "views": views,
    "max_in_memory_rows": max_rows + 1,
}

with open(str(snakemake.input["template"]), "r") as template, open(
    str(snakemake.output), "w"
) as outfile:
    result = process_yaml(template, outfile=outfile, variables=variables)
    if result is not None:
        print(result)
