from urllib.parse import urlparse
from pathlib import Path


def giga_to_byte(g):
    g = g - 2  # leave some memory to other processes
    return g * (1024 ** 3)


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
        basename = basename.replace(".gz", "")

    return basename
