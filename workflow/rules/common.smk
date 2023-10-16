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
                wildcards.serie, f"{wildcards.sample}.{ext}"
            )
            if os.path.exists(infile):
                break
    else:
        for ext in supported_extensions:
            infile = [
                trim_reads_folder.joinpath(
                    wildcards.serie, f"{wildcards.sample}_1.{ext}"
                ),
                trim_reads_folder.joinpath(
                    wildcards.serie, f"{wildcards.sample}_2.{ext}"
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
                wildcards.serie, f"{wildcards.sample}.{ext}"
            )
            if os.path.exists(candidate):
                return candidate
        raise ValueError(
            f"Could not find FastQ file. Check your naming. Supported extensions: {supported_extensions}"
        )
    elif wildcards.serie in library_names_paired:
        for ext in supported_extensions:
            candidate1 = raw_reads_folder.joinpath(
                wildcards.serie, f"{wildcards.sample}.{ext}"
            )
            candidate2 = raw_reads_folder.joinpath(
                wildcards.serie, f"{wildcards.sample}.{ext}"
            )
            if os.path.exists(candidate1) and os.path.exists(candidate2):
                return [ candidate1, candidate2 ]
        raise ValueError(
            f"Could not find FastQ files. Check your naming.\nPaired-end suffixed: {supported_suffixes}.\nSupported extensions: {supported_extensions}"
        ) 
	return ""


def get_fastq_paired(wildcards):
    if wildcards.serie in library_names_paired:
        for ext in supported_extensions:
            for suffix in supported_suffixes:
                candidate1 = raw_reads_folder.joinpath(
                    wildcards.serie, f"{wildcards.sample}{suffix[0]}.{ext}"
                )
                candidate2 = raw_reads_folder.joinpath(
                    wildcards.serie, f"{wildcards.sample}{suffix[1]}.{ext}"
                )
                if os.path.exists(candidate1) and os.path.exists(candidate2):
                    return {"m1": candidate1, "m2": candidate2}
        raise ValueError(
            f"Could not find FastQ files. Check your naming.\nPaired-end suffixed: {supported_suffixes}.\nSupported extensions: {supported_extensions}"
        )
    return ""


def mkdir(p: Path, verbose=False):
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
        if verbose:
            print("Created {}".format(p))


def get_tRNA_annotation_file(wildcards):
    checkpoint_output = checkpoints.download_gtRNAdb.get(**wildcards).output[0]
    bed_filename = glob_wildcards(os.path.join(checkpoint_output, "{x}.bed")).x[0]
    return tRNA_annotation_dir.joinpath(f"{bed_filename}.bed")


def get_trna_coverage(wildcards):
    return expand(
        trna_coverage_folder.joinpath("{{serie}}", "{sample}.bed"),
        sample=get_samples(wildcards, samples),
    )


def get_deseq2_test(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["test"]


def get_deseq2_variable(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["variable"]


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
    candidates = []
    if wildcards.serie in library_names_single:
        s = samples["single"][wildcards.serie]
        for extension in supported_extensions:
            for sample in samples["single"][wildcards.serie]:
                # we use absolute path because of Singularity/Docker
                candidates.append(raw_reads_folder.joinpath(wildcards.serie, f"{sample}.{extension}").resolve())
            if all([os.path.exists(f) for f in candidates]):
                    break
    else:
        for extension in supported_extensions:
            for m in supported_suffixes:
                candidates = []
                for sample in samples["paired"][wildcards.serie]:
                    candidates.append(raw_reads_folder.joinpath(wildcards.serie, f"{sample}{m[0]}.{extension}").resolve())
                    candidates.append(raw_reads_folder.joinpath(wildcards.serie, f"{sample}{m[1]}.{extension}").resolve())
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
            wildcards.serie, f"{wildcards.sample}.fastq.gz"
        )
    elif wildcards.serie in library_names_paired:
        return [
            trim_reads_folder.joinpath(
                wildcards.serie, f"{wildcards.sample}_1.fastq.gz"
            ),
            trim_reads_folder.joinpath(
                wildcards.serie, f"{wildcards.sample}_2.fastq.gz"
            ),
        ]
    else:
        raise ValueError(f"Could not determine protocol type for {wildcards.serie}")


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
        for m in supported_suffixes:
            for ext in supported_extensions:
                m1 = [
                    os.path.join(
                        raw_reads_folder, wildcards.serie, f"{sample}{m[0]}.{ext}"
                    )
                    for sample in s
                ]
                m2 = [
                    os.path.join(
                        raw_reads_folder, wildcards.serie, f"{sample}{m[1]}.{ext}"
                    )
                    for sample in s
                ]

                if all([os.path.exists(p) for p in m1 + m2]):
                    ret = [
                        *[
                            fastqc_raw_folder.joinpath(
                                wildcards.serie, f"{sample}{m[0]}_fastqc.html"
                            )
                            for sample in s
                        ],
                        *[
                            fastqc_raw_folder.joinpath(
                                wildcards.serie, f"{sample}{m[1]}_fastqc.html"
                            )
                            for sample in s
                        ],
                    ]

    return ret
