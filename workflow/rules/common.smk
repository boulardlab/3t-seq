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


def get_samples(wildcards, samples):
    if wildcards.serie in samples["single"]:
        s = samples["single"][wildcards.serie]
    else:
        s = samples["paired"][wildcards.serie]
    return s


def get_bw(wildcards):
    """Builds bigwig paths for rule all"""
    o = []
    for lib in library_names_single + library_names_paired:
        if lib in samples["single"].keys():
            s = samples["single"][lib]
        else:
            s = samples["paired"][lib]
        o += expand(star_folder.joinpath("{serie}", "{sample}.bw"), serie=lib, sample=s)
    return o


def get_star_input(wildcards):
    """Builds input paths for STAR alignment testing if a library is single-end or paired-end"""
    if wildcards.serie in library_names_single:
        for ext in supported_extensions:
            infile = trim_reads_folder.joinpath(
                wildcards.serie, "{0}.{1}".format(wildcards.sample, ext)
            )
            if os.path.exists(infile):
                break
    else:
        for ext in supported_extensions:
            infile = [
                trim_reads_folder.joinpath(
                    wildcards.serie, "{0}_1.{1}".format(wildcards.sample, ext)
                ),
                trim_reads_folder.joinpath(
                    wildcards.serie, "{0}_2.{1}".format(wildcards.sample, ext)
                ),
            ]
            if all([os.path.exists(f) for f in infile]):
                break
    return infile


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


def get_fastq(wildcards):
    if wildcards.serie in library_names_single:
        for ext in supported_extensions:
            candidate = raw_reads_folder.joinpath(
                wildcards.serie, "{0}.{1}".format(wildcards.sample, ext)
            )
            if os.path.exists(candidate):
                return candidate
        raise ValueError(
            "Could not find FastQ file. Check your naming. Supported extensions: {0}".format(
                supported_extensions
            )
        )
    elif wildcards.serie in library_names_paired:
        for ext in supported_extensions:
            candidate1 = raw_reads_folder.joinpath(
                wildcards.serie, "{0}.{1}".format(wildcards.sample, ext)
            )
            candidate2 = raw_reads_folder.joinpath(
                wildcards.serie, "{0}.{1}".format(wildcards.sample, ext)
            )
            if os.path.exists(candidate1) and os.path.exists(candidate2):
                return [candidate1, candidate2]
        raise ValueError(
            "Could not find FastQ files. Check your naming.\nPaired-end suffixed: {0}.\nSupported extensions: {1}".format(
                supported_suffixes, supported_extensions
            )
        )
        return ""


def get_fastq_paired(wildcards):
    if wildcards.serie in library_names_paired:
        for ext in supported_extensions:
            for suffix in supported_suffixes:
                candidate1 = raw_reads_folder.joinpath(
                    wildcards.serie,
                    "{0}{1}.{2}".format(wildcards.sample, suffix[0], ext),
                )
                candidate2 = raw_reads_folder.joinpath(
                    wildcards.serie,
                    "{0}{1}.{2}".format(wildcards.sample, suffix[1], ext),
                )
                if os.path.exists(candidate1) and os.path.exists(candidate2):
                    return {"m1": candidate1, "m2": candidate2}
        raise ValueError(
            "Could not find FastQ files. Check your naming.\nPaired-end suffixed: {0}.\nSupported extensions: {1}".format(
                supported_suffixes, supported_extensions
            )
        )
    return ""


def mkdir(p: Path, verbose=False):
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
        if verbose:
            print("Created {}".format(p))


def get_trna_coverage(wildcards):
    return expand(
        trna_coverage_folder.joinpath("{{serie}}", "{sample}.bed"),
        sample=get_samples(wildcards, samples),
    )


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
    var = deseq2_params["variable"]
    if not var in sample_sheet.columns.values.tolist():
        raise ValueError(
            "{0} was not detected in sample sheet columns. Please check your config.".format(
                var
            )
        )
    return deseq2_params["variable"]


def get_deseq2_reference_level(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["reference_level"]


def get_markdup_bam(wildcards):
    return expand(
        markdup_folder.joinpath("{{serie}}/{sample}.markdup.bam"),
        sample=get_samples(wildcards, samples),
    )


def get_markdup_fastqc(wildcards):
    return expand(
        fastqc_markdup_folder.joinpath("{{serie}}", "{sample}.markdup_fastqc.html"),
        sample=get_samples(wildcards, samples),
    )


def get_salmonTE_quant_input(wildcards):
    if wildcards.serie in library_names_single:
        s = samples["single"][wildcards.serie]
        for extension in supported_extensions:
            candidates = []
            for sample in samples["single"][wildcards.serie]:
                # we use absolute path because of Singularity/Docker
                candidates.append(
                    raw_reads_folder.joinpath(
                        wildcards.serie, "{0}.{1}".format(sample, extension)
                    ).resolve()
                )
            if all([os.path.exists(f) for f in candidates]):
                break
    else:
        for extension in supported_extensions:
            for m in supported_suffixes:
                candidates = []
                for sample in samples["paired"][wildcards.serie]:
                    candidates.append(
                        raw_reads_folder.joinpath(
                            wildcards.serie,
                            "{0}{1}.{2}".format(sample, m[0], extension),
                        ).resolve()
                    )
                    candidates.append(
                        raw_reads_folder.joinpath(
                            wildcards.serie,
                            "{0}{1}.{2}".format(sample, m[1], extension),
                        ).resolve()
                    )
                if all([os.path.exists(f) for f in candidates]):
                    break
            else:
                continue
            break
    return candidates


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


def get_star_stats(wildcards):
    return expand(
        star_folder.joinpath("{{serie}}/{sample}.Log.final.out"),
        sample=get_samples(wildcards, samples),
    )


def get_star_fastqc(wildcards):
    return expand(
        fastqc_star_folder.joinpath(
            "{{serie}}", "{sample}.Aligned.sortedByCoord.out_fastqc.html"
        ),
        sample=get_samples(wildcards, samples),
    )


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


def get_trimmomatic_stats(wildcards):
    ret = expand(
        trim_reads_folder.joinpath(wildcards.serie, "{sample}.stats.txt"),
        sample=get_samples(wildcards, samples),
    )
    return ret


def get_trimmed_fastqc(wildcards):
    s = wildcards.serie
    if s in library_names_paired:
        ret = [
            *expand(
                fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_1_fastqc.html"),
                sample=get_samples(wildcards, samples),
            ),
            *expand(
                fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_2_fastqc.html"),
                sample=get_samples(wildcards, samples),
            ),
        ]
    else:
        ret = expand(
            fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_fastqc.html"),
            sample=get_samples(wildcards, samples),
        )
    return ret


def get_fastqc(wildcards):
    s = get_samples(wildcards, samples)
    if wildcards.serie in library_names_single:
        ret = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_fastqc.html"),
            sample=s,
        )
    else:
        ret = []
        for i in range(len(supported_suffixes)):
            m = supported_suffixes[i]
            for j in range(len(supported_extensions)):
                ext = supported_extensions[j]
                m1 = []
                m2 = []
                for z in range(len(s)):
                    sample = s[z]
                    sample = sample.strip()

                    fn1 = "{0}{1}.{2}".format(sample, m[0], ext)
                    fn2 = "{0}{1}.{2}".format(sample, m[1], ext)

                    fp1 = raw_reads_folder.joinpath(wildcards.serie, fn1)
                    fp2 = raw_reads_folder.joinpath(wildcards.serie, fn2)

                    m1.append(fp1)
                    m2.append(fp2)

                if all([p.exists() for p in m1 + m2]):
                    for sample in s:
                        for i in range(len(m)):
                            ret.append(
                                fastqc_raw_folder.joinpath(
                                    wildcards.serie,
                                    "{0}{1}_fastqc.html".format(sample, m[i]),
                                )
                            )
        if not ret:
            raise ValueError("Could not determine MultiQC input files.")
    return ret


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


def get_star_counts(wildcards):
    return expand(
        star_folder.joinpath(wildcards.serie, "{sample}.ReadsPerGene.out.tab"),
        sample=get_samples(wildcards, samples),
    )
