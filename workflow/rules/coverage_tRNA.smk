
localrules:
    mk_genome_tsv,
    download_gtRNAdb,


rule mk_genome_tsv:
    input:
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        star_folder.joinpath("{serie}", "{sample}.genome"),
    conda:
        "../env/samtools.yml"
    log:
        log_folder.joinpath("mk_genome_tsv/{serie}/{sample}.log"),
    shell:
        """
        samtools view -H {input} | grep SQ | cut -f 2,3 | sed -r 's/(SN|LN)://g' | tr " " "\t" > {output}
        """


rule coverage_trna:
    input:
        bam=star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
        bai=star_folder.joinpath(
            "{serie}", "{sample}.Aligned.sortedByCoord.out.bam.bai"
        ),
        genome=star_folder.joinpath("{serie}", "{sample}.genome"),
        annotation=get_tRNA_annotation_file,
    output:
        trna_coverage_folder.joinpath("{serie}", "{sample}.bed"),
    conda:
        "../env/bedtools.yml"
    log:
        log_folder.joinpath("bedtools-trna/{serie}/{sample}.log"),
    shell:
        """
        set -x 
        SORTED=$(mktemp)
        JOINED=$(join -o 1.1 <(awk '{{print $1}}' {input.genome} | sort -u) <(awk '{{print $1}}' {input.bai} | sort -u)) 
        awk -v j="$JOINED" '$1~j' {input.genome} | grep -v 'chrUn' | bedtools sort -faidx {input.bai} -i - > $SORTED
        
        bedtools coverage -g {input.genome} -sorted -a $SORTED -b {input.bam} | sort -k1,1 -k2,2n > {output}
        rm $SORTED
        """


rule build_trna_coverage_matrix:
    input:
        get_trna_coverage,
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix.txt"),
    conda:
        "../env/R.yml"
    log:
        log_folder.joinpath("bedtools-trna/build_trna_coverage_matrix-{serie}.log"),
    script:
        "../scripts/build_tRNA_coverage_matrix_v1.R"
