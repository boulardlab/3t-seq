
localrules:
    mk_genome_tsv,
    download_gtRNAdb,


rule mk_genome_tsv:
    input:
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        star_folder.joinpath("{serie}", "{sample}.genome"),
    conda:
        "../../env/samtools.yml"
    log:
        log_folder.joinpath("mk_genome_tsv/{serie}/{sample}.log"),
    shell:
        """
        samtools view -H {input[0]} | grep SQ | cut -f 2,3 | sed -r 's/(SN|LN)://g' | tr " " "\t" > {output}
        """

def get_tRNA_annotation_file(wildcards):
    checkpoint_output = checkpoints.download_gtRNAdb.get(**wildcards).output[0]
    bed_filename = glob_wildcards(os.path.join(checkpoint_output, "{x}.bed")).x[0]
    return tRNA_annotation_dir.joinpath(f"{bed_filename}.bed")


rule coverage_trna:
    input:
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam.bai"),
        star_folder.joinpath("{serie}", "{sample}.genome"),
        get_tRNA_annotation_file,
    output:
        trna_coverage_folder.joinpath("{serie}", "{sample}.bed"),
    conda:
        "../../env/bedtools.yml"
    log:
        log_folder.joinpath("bedtools-trna/{serie}/{sample}.log"),
    shell:
        """
        set -x 
        SORTED=$(mktemp)
        JOINED=$(join -o 1.1 <(awk '{{print $1}}' {input[3]} | sort -u) <(awk '{{print $1}}' {input[2]} | sort -u)) 
        awk -v j="$JOINED" '$1~j' {input[3]} | grep -v 'chrUn' | bedtools sort -faidx {input[2]} -i - > $SORTED
        
        bedtools coverage -g {input[2]} -sorted -a $SORTED -b {input[0]} | sort -k1,1 -k2,2n > {output}
        rm $SORTED
        """


def get_trna_coverage(wildcards):
    return expand(
        trna_coverage_folder.joinpath("{{serie}}", "{sample}.bed"),
        sample=get_samples(wildcards, samples),
    )


rule build_trna_coverage_matrix:
    input:
        get_trna_coverage,
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix.txt"),
    conda:
        "../../env/R.yml"
    log:
        log_folder.joinpath("bedtools-trna/build_trna_coverage_matrix-{serie}.log"),
    script:
        "../scripts/build_tRNA_coverage_matrix_v1.R"
