import os
import json
import itertools
from yte import process_yaml

def get_name(x: str) -> str:
    bn = os.path.basename(x)
    name = os.path.splitext(bn)[0]
    return name
    

datasets = []
for dataset_path in snakemake.input["datasets"]:
    dataset_name = get_name(dataset_path)
    dataset = {
        "name": dataset_name,
        "path": dataset_path
    }
    datasets.append(dataset)

specs = [
    (spec, json.load(open(spec, "r")))
    for spec in snakemake.params["view_specs"]
]

views = []
for (dataset_path, view_spec) in itertools.product(snakemake.input["datasets"], specs):
    dataset_name = get_name(dataset_path)
    spec_name = get_name(view_spec[0])
    view_name = "-".join([dataset_name, spec_name])
    view = {
        "name": view_name,
        "dataset": dataset_name,
        "spec": view_spec[1]
    }
    views.append(view)

variables = {
    "plot_name": snakemake.params["plot_name"],
    "datasets": datasets,
    "views": views
}

with open(str(snakemake.input["template"]), "r") as template, open(
    str(snakemake.output), "w"
) as outfile:
    result = process_yaml(template, outfile=outfile, variables=variables)
    if result is not None:
        print(result)
