from pathlib import PosixPath

filepath_pattern = r"(?P<path>.*/)?(?P<sample>.+?)(?P<mate>_[MR]?[12])?(?:_sequence)?(?P<extension>\.f(?:ast)?q)(?P<gzipped>\.gz)?"
filename_pattern = r"(?P<sample>.+?)(?P<mate>_[MR]?[12])?(?:_sequence)?$"


def giga_to_byte(g):
    g = g - 2  # leave some memory to other processes
    return g * (1024**3)


def get_filename(link, decompress=False, stem=False):
    # sometimes we use this function to get filenames of paths instead of links.
    # convert back Path to str.
    if isinstance(link, Path):
        link = str(link)
    # parse
    p = urlparse(link)

    # cast to Path and get name attribute
    if stem:
        basename = Path(p.path).stem
    else:
        basename = Path(p.path).name

    # remove .gz
    if decompress and str(basename).endswith(".gz"):
        basename = str(basename).replace(".gz", "")
        basename = Path(basename)

    return basename


def get_sample_name_from_filepath(pattern, test_string):
    m = re.match(pattern, str(test_string))
    if m:
        gd = m.groupdict()
        return gd["sample"]
    else:
        raise ValueError(
            "File path does not match expected pattern: {}".format(test_string)
        )


def get_sample_sheet_path(wildcards):
    sample_sheet_path = next(
        (
            lib["sample_sheet"]
            for lib in config["sequencing_libraries"]
            if lib["name"] == wildcards.serie
        ),
        None,
    )
    sample_sheet_path = Path(sample_sheet_path)
    return sample_sheet_path


def get_samples(wildcards) -> list[str]:

    sample_sheet_path = get_sample_sheet_path(wildcards)

    samples = list()

    sample_sheet = pd.read_csv(sample_sheet_path)

    colname = "filename"
    if wildcards.serie in library_names_paired:
        colname += "_1"

    for fn in sample_sheet[colname].tolist():
        sample_name = get_sample_name_from_filepath(filename_pattern, str(fn))
        samples.append(sample_name)

    return samples


def get_samples_names(wildcards) -> list[str]:
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)
    return sample_sheet["name"].tolist()


def get_bw(wildcards):
    """Builds bigwig paths for rule all"""
    o = []
    for lib in config["sequencing_libraries"]:
        sample_sheet = pd.read_csv(lib["sample_sheet"])
        samples = sample_sheet["name"].tolist()

        o += expand(
            star_folder.joinpath("{serie}", "{sample}.bw"),
            serie=lib["name"],
            sample=list(samples),
        )
    return o


def get_star_input(wildcards):
    """Builds input paths for STAR alignment testing if a library is single-end or paired-end"""

    samples_names = get_samples_names(wildcards)
    filenames_no_mate = get_samples(wildcards)
    idx = samples_names.index(wildcards.sample)
    s = filenames_no_mate[idx]

    if wildcards.serie in library_names_single:
        ret = trim_reads_folder.joinpath(wildcards.serie, "{0}.fastq.gz".format(s))
    else:
        ret = [
            trim_reads_folder.joinpath(wildcards.serie, "{0}_1.fastq.gz".format(s)),
            trim_reads_folder.joinpath(wildcards.serie, "{0}_2.fastq.gz".format(s)),
        ]
    return ret


def get_params(wildcards, key):
    """Returns the value of a specific key for the current serie"""
    params = ""
    for lib in config["sequencing_libraries"]:
        if lib["name"] == wildcards.serie:
            params = lib[key]
    return params


def get_sample_sheet(wildcards):
    """Returns path to sample sheet for current serie"""
    return get_params(wildcards, "sample_sheet")


def parse_filepath(filepath: PosixPath):
    m = re.match(filepath_pattern, str(filepath))
    if m:
        gd = m.groupdict()
        return gd
    raise ValueError("Invalid path: {}".format(filepath))


def get_fastq_paired(wildcards):

    ret = {"m1": "", "m2": ""}

    for p in raw_reads_folder.joinpath(wildcards.serie).iterdir():
        gd = parse_filepath(p)
        sample = gd["sample"]
        mate = gd["mate"]

        if re.search(sample, wildcards.sample):
            if re.search(r"1", mate):
                ret["m1"] = p
            elif re.search(r"2", mate):
                ret["m2"] = p
            else:
                raise ValueError(
                    "Could not find paired input files for sample {}.\nMate: {}\nFull path: {}".format(
                        wildcards.sample, mate, str(p)
                    )
                )
    return ret


def get_fastq(wildcards):
    for p in raw_reads_folder.joinpath(wildcards.serie).iterdir():
        gd = parse_filepath(p)
        s = gd["sample"]
        if gd["mate"] != "":
            s += gd["mate"]
        if s == wildcards.sample:
            return p
    raise ValueError(
        "Could not determine input files for serie: {}\nsample: {}\n{}".format(
            wildcards.serie, wildcards.sample, gd
        )
    )


def mkdir(p: Path, verbose=False):
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
        if verbose:
            print("Created {}".format(p))


def get_trna_coverage(wildcards):
    sample_sheet_path = get_sample_sheet_path(wildcards)
    return {
        "bed": expand(
            trna_coverage_folder.joinpath(wildcards.serie, "{sample}.bed"),
            sample=get_samples_names(wildcards),
        ),
        "sample_sheet": sample_sheet_path,
    }


def get_deseq2_test(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    if not "test" in deseq2_params:
        deseq2_params["test"] = "Wald"
    test = deseq2_params["test"]
    if not test in ["Wald", "LRT"]:
        raise ValueError(
            "Invalid test: {0}. Test name must be either Wald or LRT. Check your config.".format(
                test
            )
        )
    return deseq2_params["test"]


def get_deseq2_variable(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    if not "variable" in deseq2_params:
        deseq2_params["variable"] = "genotype"
    # var = deseq2_params["variable"]
    # if not var in sample_sheet.columns.values.tolist():
    #     raise ValueError(
    #         "{0} was not detected in sample sheet columns. Please check your config.".format(
    #             var
    #         )
    #     )
    return deseq2_params["variable"]


def get_deseq2_reference_level(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["reference_level"]


def get_markdup_bam(wildcards):
    return {
        "bam": expand(
            markdup_folder.joinpath("{{serie}}/{sample}.markdup.bam"),
            sample=get_samples_names(wildcards),
        ),
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def get_markdup_fastqc(wildcards):
    return {
        "fastqc": expand(
            fastqc_markdup_folder.joinpath("{{serie}}", "{sample}.markdup_fastqc.zip"),
            sample=get_samples_names(wildcards),
        ),
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def get_salmonTE_quant_input(wildcards):
    ret = []
    for p in raw_reads_folder.joinpath(wildcards.serie).iterdir():
        if re.match(filepath_pattern, str(p)):
            ret.append(p.resolve())
    return ret


def set_salmonTE_genome():
    genome_label = config["genome"]["label"][:2]
    salmon_label = ""
    if genome_label == "mm":
        salmon_label = "mm"
    elif genome_label == "hg":
        salmon_label = "hs"
    elif genome_label == "dr":
        salmon_label = "dr"
    elif genome_label == "dm":
        salmon_label = "dm"
    else:
        raise ValueError(f'Unsupported genome label: {config["genome"]["label"]}')
    return salmon_label


def get_multiqc_star_inputs(wildcards):
    return {
        "star_stats": expand(
            star_folder.joinpath("{{serie}}/{sample}.Log.final.out"),
            sample=get_samples_names(wildcards),
        ),
        "fastqc": expand(
            fastqc_star_folder.joinpath(
                "{{serie}}", "{sample}.Aligned.sortedByCoord.out_fastqc.zip"
            ),
            sample=get_samples_names(wildcards),
        ),
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def get_trimmed_fastq(wildcards):
    if wildcards.serie in library_names_single:
        return trim_reads_folder.joinpath(
            wildcards.serie, "{0}.fastq.gz".format(wildcards.sample)
        )
    elif wildcards.serie in library_names_paired:
        return [
            trim_reads_folder.joinpath(
                wildcards.serie, "{0}_1.fastq.gz".format(wildcards.sample)
            ),
            trim_reads_folder.joinpath(
                wildcards.serie, "{1}_2.fastq.gz".format(wildcards.sample)
            ),
        ]
    else:
        raise ValueError(
            "Could not determine protocol type for {0}".format(wildcards.serie)
        )


def get_multiqc_trim_inputs(wildcards):
    s = wildcards.serie
    if s in library_names_paired:
        fastqcs = [
            *expand(
                fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_1_fastqc.zip"),
                sample=get_samples(wildcards),
            ),
            *expand(
                fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_2_fastqc.zip"),
                sample=get_samples(wildcards),
            ),
        ]
    else:
        fastqcs = expand(
            fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_fastqc.zip"),
            sample=get_samples(wildcards),
        )

    return {
        "trimmomatic_stats": expand(
            trim_reads_folder.joinpath(wildcards.serie, "{sample}.stats.txt"),
            sample=get_samples(wildcards),
        ),
        "fastqc": fastqcs,
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def get_fastqc(wildcards):
    s = get_samples(wildcards)
    if wildcards.serie in library_names_single:
        fastqcs = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_fastqc.zip"),
            sample=s,
        )
    else:
        mates = set()
        for p in raw_reads_folder.joinpath(wildcards.serie).iterdir():
            gd = parse_filepath(p)
            mates.add(gd["mate"])
        fastqcs = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}{mate}_fastqc.zip"),
            sample=s,
            mate=mates,
        )
    return {"fastqc": fastqcs, "sample_sheet": get_sample_sheet_path(wildcards)}


def build_rule_all_inputs(wildcards):
    ret = []

    # MultiQC reports at different steps
    ret += expand(
        multiqc_raw_folder.joinpath("{serie}", "multiqc_report.html"),
        serie=library_names_single + library_names_paired,
    )

    ret += expand(
        multiqc_trim_folder.joinpath("{serie}", "multiqc_report.html"),
        serie=library_names_single + library_names_paired,
    )
    ret += expand(
        multiqc_star_folder.joinpath("{serie}", "multiqc_report.html"),
        serie=library_names_single + library_names_paired,
    )
    ret += expand(
        multiqc_markdup_folder.joinpath("{serie}", "multiqc_report.html"),
        serie=library_names_single + library_names_paired,
    )

    # DESeq2 flags
    ret += expand(
        rdata_folder.joinpath("deseq2/{serie}/dds.rds"),
        serie=library_names_single + library_names_paired,
    )
    ret += expand(
        analysis_folder.joinpath("datavzrd", "{serie}", "datavzrd"),
        serie=library_names_single + library_names_paired,
    )
    # Bigwig files
    ret += get_bw(wildcards)

    # TE
    if not config["disable_TE_analysis"]:
        # SalmonTE results folders
        ret += expand(
            data_folder.joinpath("salmonTE/de_analysis/{se_serie}"),
            se_serie=library_names_single,
        )
        ret += expand(
            data_folder.joinpath("salmonTE/de_analysis/{pe_serie}"),
            pe_serie=library_names_paired,
        )
        # FeatureCounts tables from STAR-TE
        ret += expand(
            starTE_folder.joinpath("{se_serie}/featureCount/{method}.txt"),
            se_serie=library_names_single,
            method=["multihit", "random"],
        )
        ret += expand(
            starTE_folder.joinpath("{pe_serie}/featureCount/{method}.txt"),
            pe_serie=library_names_paired,
            method=["multihit", "random"],
        )
        ret += expand(
            starTE_folder.joinpath("{serie}", "random", "datavzrd"),
            serie=library_names_single + library_names_paired,
        )

    # tRNA coverage files
    if not config["disable_tRNA_analysis"]:
        ret += expand(
            trna_coverage_folder.joinpath("{serie}", "tRNA_lfc.txt"),
            serie=library_names_paired + library_names_single,
        )
        ret += expand(
            trna_coverage_folder.joinpath("{serie}", "datavzrd"),
            serie=library_names_paired + library_names_single,
        )

    return ret


def get_deseq2_inputs(wildcards):
    return {
        "star_counts": expand(
            star_folder.joinpath(wildcards.serie, "{sample}.ReadsPerGene.out.tab"),
            sample=get_samples_names(wildcards),
        ),
        "annotation_file": gtf_path,
        "sample_sheet": get_sample_sheet_path(wildcards),
    }
