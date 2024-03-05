fasta_file = open(snakemake.input["genome_fasta_file"], "r")
annotation_file = open(snakemake.input["genome_annotation_file"], "r")


def test_file(fh):
    ret = True
    while i < 1000:
        line = fh.readline()
        if not line.startswith(">") and not line.startswith(">chr"):
            ret = False
    fh.seek(0, 0)
    return ret


fasta_has_chr = test_file(fasta_file)
annotation_has_chr = test_file(annotation_file)
