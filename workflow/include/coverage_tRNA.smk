from utilities.rnaseq import get_samples

localrules: mk_genome_tsv, download_gtRNAdb

rule download_gtRNAdb:
    output:
        tRNA_annotation_file
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    shell:
        """
        wget -O http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz
        tar xf mm10-tRNAs.tar.gz
        mv mm10-tRNAs.bed {output}
        rm mm10-*
        """


rule mk_genome_tsv:
    input:
        star_folder.joinpath("{serie}","{sample}.Aligned.sortedByCoord.out.bam")
    output:
        star_folder.joinpath("{serie}","{sample}.genome")
    singularity:
        str(container_folder.joinpath("samtools.sif"))
    log:
        log_folder.joinpath("mk_genome_tsv/{serie}/{sample}.log")
    shell:
        """
        samtools view -H {input[0]} | grep SQ | cut -f 2,3 | sed -r 's/(SN|LN)://g' | tr " " "\t" > {output}
        """

rule coverage_trna:
    input:
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
        star_folder.joinpath("{serie}","{sample}.Aligned.sortedByCoord.out.bam.bai"),
        star_folder.joinpath("{serie}","{sample}.genome"),
        tRNA_annotation_file
    output:
        trna_coverage_folder.joinpath("{serie}", "{sample}.bed")
    singularity:
        str(container_folder.joinpath("bedtools.sif"))
    log:
        log_folder.joinpath("bedtools-trna/{serie}/{sample}.log")
    shell:
        """
        SORTED=$(mktemp)
        grep -v 'chrUn' {input[3]} | bedtools sort -faidx {input[2]} -i - > $SORTED        
        bedtools coverage -g {input[2]} -sorted -a $SORTED -b {input[0]} | sort -k1,1 -k2,2n > {output}
        rm $SORTED
        """

def get_trna_coverage(wildcards):
    return expand(trna_coverage_folder.joinpath("{{serie}}", "{sample}.bed"), sample = get_samples(wildcards, samples))

rule build_trna_coverage_matrix:
    input:
        get_trna_coverage
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix.txt")
    singularity:
        str(container_folder.joinpath("R.sif"))
    log:
        log_folder.joinpath("bedtools-trna/build_trna_coverage_matrix-{serie}.log")
    script:
        "../../../src/R/build_tRNA_coverage_matrix_v1.R"

