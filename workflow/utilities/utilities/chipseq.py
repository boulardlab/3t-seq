
def build_sample_wildcard(wildcards, samples, wildcard_name, index):
    if wildcards.serie in samples["pe"]:
        case_control = samples["pe"][wildcards.serie]["case-control"]
    else:
        case_control = samples["se"][wildcards.serie]["case-control"]
    s = [x[index] for x in case_control]
    s = list(set(s))
    w = "|".join(s)
    w = "{"+ wildcard_name +"," + w + "}"
    return w

