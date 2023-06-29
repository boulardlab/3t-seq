
def get_samples(wildcards, samples):
    if wildcards.serie in samples["single"]:
        s = samples["single"][wildcards.serie]
    else:
        s = samples["paired"][wildcards.serie]
    return s
