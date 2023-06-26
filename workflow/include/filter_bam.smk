rule filter_bam:
    input:
        alignment=starTE_folder.joinpath("{serie}/{method}/{sample}.Aligned.out.bam"),
        annotation=rmsk_bed
    output:
        starTE_folder.joinpath("{serie}/filter/{method}/{sample}.TEonly.bam")
    log:
        log_folder.joinpath("samtools_view/{serie}/{method}/{sample}.log")
    threads: 2
    singularity:
        str(container_folder.joinpath("samtools.sif"))
    shell:
        """
        set -x
        samtools view -L {input.annotation} -@ {threads} -o {output} {input.alignment} |& tee {log}
        """

# rule filter_bam_random:
#     input:
#         alignment=starTE_folder.joinpath("{serie}/random/{sample}.Aligned.out.bam"),
#         annotation=rmsk_bed
#     output:
#         starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam")
#     log:
#         log_folder.joinpath("samtools_view/{serie}/random/{sample}.log")
#     threads: 2
#     singularity:
#         str(container_folder.joinpath("samtools.sif"))
#     shell:
#         """
#         set -x
#         samtools view -L {input.annotation} -@ {threads} {input.alignment} |& tee {log}
#         """
