rule featureCounts_multihit:
    input:
        bam=lambda wildcards: expand(starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam"),
                                     serie=wildcards.serie, sample=samples['single'][wildcards.serie] if wildcards.serie in samples['single'] else samples['paired'][wildcards.serie]),
        annotation=rmsk_path
    output:
        starTE_folder.joinpath("{serie}/featureCount/multihit.txt")
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
        log_folder.joinpath("featureCounts/{serie}/multihit.log")
    threads:
        4
    shell:
         """
         set -x
         featureCounts -M --fraction -F GTF -T {threads} -s 0 -a {input.annotation} -o {output} {input.bam}
         """
