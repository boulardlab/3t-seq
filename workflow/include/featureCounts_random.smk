rule featureCounts_random:
    input:
        bam=lambda wildcards: expand(starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam"),
                                     serie=wildcards.serie, sample=samples['single'][wildcards.serie] if wildcards.serie in samples['single'] else samples['paired'][wildcards.serie]),
        annotation=rmsk_path
    output:
        starTE_folder.joinpath("{serie}/featureCount/random.txt")
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
        log_folder.joinpath("featureCounts/{serie}/random.log")
    threads:
        4
    shell:
         """
         set -x
         featureCounts -M -F GTF -T {threads} -s 0 -a {input.annotation} -o {output} {input.bam}
         """
