from urllib import request, parse
import json
import csv
import sys


def main(
    genome="mm10", track_name="rmsk", gtf_output=sys.stdout, bed_output=sys.stdout
):
    """
    A function that queries UCSC genome browser to download RepeatMasker
    annotation.
    """

    # setup query string
    query_object = {"genome": genome, "track": track_name}
    query_string = parse.urlencode(query_object).encode("ascii")

    # output field names
    gtf_fieldnames = [
        "genoName",
        "source",
        "feature",
        "genoStart",
        "genoEnd",
        "swScore",
        "strand",
        "frame",
        "attribute",
    ]

    bed_fieldnames = ["genoName", "genoStart", "genoEnd", "swScore", "strand", "name"]

    # create request object
    url = request.Request(
        url="https://api.genome.ucsc.edu/getData/track",
        data=query_string,
        method="GET",
    )

    # create DictWriter for later
    gtf_writer = csv.DictWriter(
        f=gtf_output,
        fieldnames=gtf_fieldnames,
        dialect="excel-tab",
    )

    bed_writer = csv.DictWriter(
        f=bed_output, fieldnames=bed_fieldnames, dialect="excel-tab"
    )

    # perform the query and parse response
    with request.urlopen(url) as response:
        dat = json.loads(response.read())
        track_data = dat["rmsk"]

        for chromosome in track_data.values():
            for repeat in chromosome:
                # UCSC works zero-based base position
                # https://www.biostars.org/p/84686/

                # GTF is 1-based, therefore, genoStart needs +1
                gtf_row = {
                    "genoName": repeat["genoName"],
                    "source": "RepeatMasker",
                    "feature": "exon",
                    "genoStart": int(repeat["genoStart"]) + 1,
                    "genoEnd": int(repeat["genoEnd"]),
                    "swScore": repeat["swScore"],
                    "strand": repeat["strand"],
                    "frame": ".",
                    "attribute": ";".join(
                        [
                            f"{key}={repeat[key]}"
                            for key in filter(
                                lambda x: x
                                in [
                                    "repName",
                                    "repClass",
                                    "repFamily",
                                    "repStart",
                                    "repEnd",
                                    "repLeft",
                                ],
                                repeat,
                            )
                        ]
                    ),
                }

                # BED is zero-based, therefore, genoStart does not need +1
                bed_row = {
                    "genoName": repeat["genoName"],
                    "genoStart": repeat["genoStart"],
                    "genoEnd": repeat["genoEnd"],
                    "swScore": repeat["swScore"],
                    "strand": repeat["strand"],
                    "name": repeat["repName"],
                }

                gtf_writer.writerow(gtf_row)
                bed_writer.writerow(bed_row)


if __name__ == "__main__":
    genome = snakemake.params["genome_id"]
    track_name = "rmsk"
    gtf_output = open(str(snakemake.output[0]), "w")
    bed_output = open(str(snakemake.output[1]), "w")

    main(genome, track_name, gtf_output, bed_output)
