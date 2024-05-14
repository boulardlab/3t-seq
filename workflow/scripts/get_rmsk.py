from urllib import request, parse
from urllib.error import HTTPError, URLError
import json
import csv
import sys
import logging
from time import sleep
from random import randint
import re


def setup_logger(log_file):
    # Create a logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # Create a file handler and set the level to debug
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)

    # Create a console handler and set the level to info
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Create a formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def build_request(base_url, params, logger):
    query_string = parse.urlencode(params).encode("ascii")

    url = request.Request(
        url=base_url,
        data=query_string,
        method="GET",
    )
    if logger:
        logger.debug("Created url: {0}".format(url.get_full_url()))
    return url


def fetch_chromosomes(genome, track, logger):
    """
    A function to fetch the list of chromosomes of a given organism from UCSC.
    """

    base_url = "https://api.genome.ucsc.edu/list/chromosomes"

    params = {"genome": genome, "track": track}
    url = build_request(base_url, params, logger)

    try:
        with request.urlopen(url) as response:
            data = json.loads(response.read().decode("utf-8"))
            data["chromosomes"] = list(
                filter(
                    lambda chrom: re.match(r"^chr(?:[1-1][0-9]|X|Y|M)$", chrom),
                    data["chromosomes"],
                )
            )
            return data
    except HTTPError as e:
        if logger:
            logger.critical(f"HTTPError: {e.code} - {e.reason}")
        raise e
    except URLError as e:
        if logger:
            logger.critical(f"URLError: {e.reason}")
        raise e


def get_gtf_writer(gtf_output, logger):
    """
    A function to build a GTF-compliant CSV writer.
    """

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

    # create DictWriter for later
    gtf_writer = csv.DictWriter(
        f=gtf_output,
        fieldnames=gtf_fieldnames,
        dialect="excel-tab",
    )
    if logger:
        logger.debug("Created GTF writer. Fields = {0}".format(gtf_fieldnames))

    return gtf_writer


def get_gtf_row(data):
    """
    A function that maps UCSC RepeatMasker repeat JSON to GTF fields
    """
    # GTF is 1-based, therefore, genoStart needs +1
    return {
        "genoName": data["genoName"],
        "source": "RepeatMasker",
        "feature": "exon",
        "genoStart": int(data["genoStart"]) + 1,
        "genoEnd": int(data["genoEnd"]),
        "swScore": data["swScore"],
        "strand": data["strand"],
        "frame": ".",
        "attribute": ";".join(
            [
                f"{key}={data[key]}"
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
                    data,
                )
            ]
        ),
    }


def write_gtf_row(writer, data):
    """
    A function that writes a single GTF row using the given GTF Writer and data
    """
    row = get_gtf_row(data)
    writer.writerow(row)


def get_bed_writer(bed_output, logger):
    """
    A function to build a BED-compliant CSV writer.
    """
    bed_fieldnames = ["genoName", "genoStart", "genoEnd", "swScore", "strand", "name"]
    bed_writer = csv.DictWriter(
        f=bed_output, fieldnames=bed_fieldnames, dialect="excel-tab"
    )
    if logger:
        logger.debug("Created BED writer. Fields = {0}".format(bed_fieldnames))
    return bed_writer


def get_bed_row(data):
    """
    A function that maps UCSC RepeatMasker repeat JSON to BED fields
    """
    # BED is zero-based, therefore, genoStart does not need +1
    return {
        "genoName": data["genoName"],
        "genoStart": data["genoStart"],
        "genoEnd": data["genoEnd"],
        "swScore": data["swScore"],
        "strand": data["strand"],
        "name": data["repName"],
    }


def write_bed_row(writer, data):
    """
    A function that writes a single BED row using the given BED Writer and data
    """
    row = get_bed_row(data)
    writer.writerow(row)


def main(
    genome="mm10",
    track_name="rmsk",
    gtf_output=sys.stdout,
    bed_output=sys.stdout,
    selected_chromosomes=None,
    logger=None,
):
    """
    A function that queries UCSC genome browser to download RepeatMasker
    annotation.
    """

    gtf_writer = get_gtf_writer(gtf_output, logger)
    bed_writer = get_bed_writer(bed_output, logger)

    if not selected_chromosomes:
        chromosome_dict = fetch_chromosomes(genome, track_name, logger)
    else:
        chromosome_dict = {"chromosomes": selected_chromosomes}

    base_url = "https://api.genome.ucsc.edu/getData/track"

    counter = 0
    nchromosomes = 0

    for chrom in chromosome_dict["chromosomes"]:
        nchromosomes += 1

        # setup query string
        params = {"genome": genome, "track": track_name, "chrom": chrom}
        if logger:
            logger.debug("Query string = {0}".format(params))

        # create request object
        url = build_request(base_url, params, logger)

        if logger:
            logger.debug("Starting GTF/BED export for chromosome {0}".format(chrom))

        chrom_lines = 0
        try:
            # perform the query and parse response
            with request.urlopen(url) as response:
                dat = json.loads(response.read())

                if logger:
                    logger.info(
                        "Downloaded {0} - {1} repeats - track name {2} - url {3}".format(
                            chrom, len(dat["rmsk"]), track_name, url.get_full_url()
                        )
                    )

                for repeat in dat["rmsk"]:
                    # UCSC works zero-based base position
                    # https://www.biostars.org/p/84686/

                    write_gtf_row(gtf_writer, repeat)
                    write_bed_row(bed_writer, repeat)

                    chrom_lines += 1
                    counter += 1

            if logger:
                logger.debug(
                    "Completed {0} - written {1} lines.".format(chrom, chrom_lines)
                )

            sleep(randint(1, 3))

        except HTTPError as e:
            if logger:
                logger.critical(f"HTTPError: {e.code} - {e.reason}")
            raise e
        except URLError as e:
            if logger:
                logger.critical(f"URLError: {e.reason}")
            raise e

    if logger:
        logger.info(
            "Successfully written {0} lines over {1} chromosomes.".format(
                counter, nchromosomes
            )
        )


if __name__ == "__main__":
    logger = setup_logger(str(snakemake.log))

    genome = snakemake.params["genome_id"]
    logger.info("Genome = {0}".format(genome))

    track_name = "rmsk"
    logger.info("Track = {0}".format(track_name))

    gtf_output = open(str(snakemake.output[0]), "w")
    logger.info("Output path = {0}".format(gtf_output))

    bed_output = open(str(snakemake.output[1]), "w")
    logger.info("Bed output = {0}".format(bed_output))

    selected_chromosomes = None
    if snakemake.params["selected_chromosome"] and isinstance(
        snakemake.params["selected_chromosome"], list
    ):
        logger.debug('Detected snakemake.params["selected_chromosome"] is a list.')

        selected_chromosomes = snakemake.params["selected_chromosome"]
        logger.info("Selected {0} for downloading".format(selected_chromosomes))
    elif snakemake.params["selected_chromosome"] and isinstance(
        snakemake.params["selected_chromosome"], str
    ):
        logger.debug(
            'Detected snakemake.params["selected_chromosome"] is a string. Converting to list.'
        )

        selected_chromosomes = [snakemake.params["selected_chromosome"]]
        logger.info("Selected {0} for downloading".format(selected_chromosomes))
    else:
        selected_chromosomes = None
        logger.info("Downloading all chromosomes")

    main(
        genome=genome,
        track_name=track_name,
        gtf_output=gtf_output,
        bed_output=bed_output,
        selected_chromosomes=selected_chromosomes,
        logger=logger,
    )
