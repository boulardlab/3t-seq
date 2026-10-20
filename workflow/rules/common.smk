import hashlib
import json
import yaml
import copy

from pathlib import Path, PosixPath
from jsonschema import Draft7Validator, validators
from snakemake.io import Wildcards
from snakemake.exceptions import WorkflowError


class Struct:
    """A simple class to convert a dictionary into an object with attributes."""

    def __init__(self, **entries):
        self.__dict__.update(entries)


filepath_pattern = r"(?P<path>.*/)?(?P<sample>.+?)(?P<mate>_[MRr]?[12])?(?P<genecore_suffix>_sequence)?(?P<extension>\.f(?:ast)?q)(?P<gzipped>\.gz)?"
filename_pattern = r"(?P<sample>.+?)(?P<mate>_[MRr]?[12])?(?:_sequence)?(?P<extension>\.f(?:ast)?q)?(?P<gzipped>\.gz)?$"


def apply_schema_defaults(config, schema_path):
    """
    Parses a JSON schema file and populates the config dictionary with
    missing values that have a 'default' defined in the schema.
    Uses the official jsonschema extension pattern to handle nested defaults correctly.
    """

    with open(schema_path, "r") as f:
        schema = yaml.safe_load(f)

    def extend_with_default(validator_class):
        validate_properties = validator_class.VALIDATORS["properties"]

        def set_defaults(validator, properties, instance, schema):
            for property, subschema in properties.items():
                if "default" in subschema:
                    if property not in instance or instance[property] is None:
                        instance[property] = copy.deepcopy(subschema["default"])

            for error in validate_properties(validator, properties, instance, schema):
                yield error

        return validators.extend(
            validator_class,
            {"properties": set_defaults},
        )

    # Use specialized validator to populate defaults recursively
    DefaultValidatingDraft7Validator = extend_with_default(Draft7Validator)
    DefaultValidatingDraft7Validator(schema).validate(config)


MIGRATION_GUIDE = (
    "docs/migration/ (https://github.com/boulardlab/3t-seq/tree/main/docs/migration)"
)

# Genomes whose references (refgenie FASTA/GTF, UCSC RepeatMasker, GtRNAdb)
# are fetched automatically. Keep in sync with the genome.label anyOf in the schema.
AUTO_FETCH_GENOMES = ("mm10", "mm39")

_DERIVED_FOLDER_HINT = (
    "Output folders are now derived from 'globals.results_folder'; remove this key."
)

# Config keys removed or moved since v1.1.0, keyed by path ("*" matches a list
# index). The schema rejects them (additionalProperties: false); this table
# only turns the generic error into an actionable one.
REMOVED_CONFIG_KEYS = {
    ("deseq2",): (
        "The top-level 'deseq2' block was moved: set test/variable/reference_level "
        "under 'defaults.deseq2' (or per comparison under 'comparisons[].deseq2'). "
        "'working_directory' and 'notebook_path' were removed."
    ),
    ("globals", "reads_folder"): (
        "Raw reads are now located through the sample sheet: put paths (absolute, "
        "or relative to the working directory) in the 'filename' (se) or "
        "'filename_1' and 'filename_2' (pe) columns."
    ),
    ("globals", "qc_folder"): _DERIVED_FOLDER_HINT,
    ("globals", "log_folder"): _DERIVED_FOLDER_HINT,
    ("globals", "references_folder"): _DERIVED_FOLDER_HINT,
    ("globals", "analysis_folder"): _DERIVED_FOLDER_HINT,
    ("globals", "data_folder"): _DERIVED_FOLDER_HINT,
    ("globals", "container_folder"): (
        "Containers are configured through the Snakemake profile; remove this key."
    ),
    ("globals", "tmp_folder"): (
        "Temporary files go under 'globals.results_folder' or $TMPDIR; remove this key."
    ),
    ("genome", "fasta_url"): (
        "Remote FASTA URLs are no longer supported: download the file and set "
        f"'genome.fasta_path', or use label {' or '.join(AUTO_FETCH_GENOMES)} to fetch it automatically."
    ),
    ("genome", "fasta_link"): (
        "Remote FASTA URLs are no longer supported: download the file and set "
        f"'genome.fasta_path', or use label {' or '.join(AUTO_FETCH_GENOMES)} to fetch it automatically."
    ),
    ("genome", "gtf_url"): (
        "Remote GTF URLs are no longer supported: download the file and set "
        f"'genome.gtf_path', or use label {' or '.join(AUTO_FETCH_GENOMES)} to fetch it automatically."
    ),
    ("genome", "gaf_url"): "No longer used; remove this key.",
    ("genome", "rmsk_path"): (
        "Custom RepeatMasker annotations are no longer supported; it is downloaded "
        f"automatically for {' or '.join(AUTO_FETCH_GENOMES)}."
    ),
    ("genome", "gtrnadb_url"): (
        "Custom GtRNAdb annotations are no longer supported; they are downloaded "
        f"automatically for {' or '.join(AUTO_FETCH_GENOMES)}."
    ),
    ("genome", "gtrnadb_bed"): (
        "Custom GtRNAdb annotations are no longer supported; they are downloaded "
        f"automatically for {' or '.join(AUTO_FETCH_GENOMES)}."
    ),
}

# Config keys renamed since v1.1.0 that are still accepted: old -> new.
DEPRECATED_CONFIG_KEYS = {"sequencing_libraries": "comparisons"}


def _format_config_path(path) -> str:
    """Render a jsonschema error path as e.g. 'comparisons[0].protocol'."""
    out = ""
    for part in path:
        out += f"[{part}]" if isinstance(part, int) else (f".{part}" if out else part)
    return out or "<config root>"


def _describe_config_error(error, config: dict) -> str:
    """Turn a jsonschema ValidationError into an actionable message.

    Args:
        error: a jsonschema.ValidationError raised against `config`.
        config: the validated config, used to name comparisons in messages.

    Returns:
        A single (possibly multi-line) message string.
    """
    path = list(error.absolute_path)
    where = _format_config_path(path)
    if path[:1] == ["comparisons"] and len(path) > 1 and isinstance(path[1], int):
        entry = config["comparisons"][path[1]]
        name = entry.get("name") if isinstance(entry, dict) else None
        if name:
            where += f" (comparison '{name}')"

    if error.validator == "additionalProperties":
        allowed = error.schema.get("properties", {})
        extra = [k for k in error.instance if k not in allowed]
        lines = []
        for key in extra:
            normalized = tuple("*" if isinstance(p, int) else p for p in path + [key])
            hint = REMOVED_CONFIG_KEYS.get(
                normalized,
                f"Unknown key. Allowed keys here: {', '.join(sorted(allowed))}.",
            )
            lines.append(f"{_format_config_path(path+[key])}: {hint}")
        return "\n  ".join(lines)

    if error.validator == "required" and path[:1] == ["comparisons"]:
        if "'protocol'" in error.message:
            return (
                f"{where}: missing 'protocol'. Set 'protocol: pe' (paired-end, sample "
                "sheet columns filename_1/filename_2) or 'protocol: se' (single-end, "
                "column filename). 1.x inferred it from the sample sheet; 2.0 does not."
            )

    if error.validator == "anyOf" and path == ["genome"]:
        return (
            f"genome.label: '{error.instance.get('label')}' has no automatically "
            f"fetched references (supported: {', '.join(AUTO_FETCH_GENOMES)}). "
            "For any other genome set 'genome.fasta_path', 'genome.gtf_path' and "
            "'genome.annotation_type'."
        )

    if error.validator == "pattern" and path == ["genome", "label"]:
        return (
            f"genome.label: '{error.instance}' is not a safe label. Use 1-64 letters, "
            "digits, '.', '_' or '-', starting with a letter or digit (it is used in file names)."
        )

    if error.validator == "dependencies" and path == ["genome"]:
        return (
            f"genome: {error.message}. Custom references need 'fasta_path', "
            "'gtf_path' and 'annotation_type' together."
        )

    return f"{where}: {error.message}"


def check_genome_support(config: dict) -> None:
    """Fail early when enabled analyses need references not available for the genome.

    Only mm10/mm39 have RepeatMasker and GtRNAdb annotations fetched automatically.
    SalmonTE is handled by disable_unsupported_salmonTE() instead.

    Raises:
        WorkflowError: listing each enabled analysis that cannot run and how to
            disable it.
    """
    label = config["genome"]["label"]
    if label in AUTO_FETCH_GENOMES:
        return

    problems = []
    if not config["disable_TE_analysis"]:
        problems.append(
            "STAR-TE analysis needs a RepeatMasker annotation, which is only fetched "
            f"for {' or '.join(AUTO_FETCH_GENOMES)}. Set 'disable_TE_analysis: true'."
        )
    if not config["disable_tRNA_analysis"]:
        problems.append(
            "tRNA analysis needs GtRNAdb annotations, which are only fetched for "
            f"{' or '.join(AUTO_FETCH_GENOMES)}. Set 'disable_tRNA_analysis: true'."
        )

    if problems:
        raise WorkflowError(
            f"genome.label '{label}' is not one of {', '.join(AUTO_FETCH_GENOMES)}, "
            "so these enabled analyses cannot run:\n  - "
            + "\n  - ".join(problems)
            + "\nIf your custom references are an mm10/mm39 assembly, use that label instead."
        )


def get_salmonTE_genome(label: str) -> str | None:
    """Return the SalmonTE reference genome for a genome label, or None if unsupported."""
    return SALMONTE_GENOMES.get(label[:2])


def disable_unsupported_salmonTE(config: dict) -> None:
    """Turn SalmonTE off, with one warning, when it has no index for the genome.

    An explicit 'disable_salmonTE_analysis: true' stays silent.
    """
    if config["disable_salmonTE_analysis"]:
        return
    label = config["genome"]["label"]
    if get_salmonTE_genome(label) is not None:
        return

    config["disable_salmonTE_analysis"] = True
    if workflow.is_main_process:
        logger.warning(
            f"SalmonTE DISABLED: genome.label '{label}' does not start with a prefix "
            f"SalmonTE has a TE index for ({', '.join(SALMONTE_GENOMES)}), so all "
            "SalmonTE rules and outputs are skipped. Set "
            "'disable_salmonTE_analysis: true' to silence this warning."
        )


def validate_config(config: dict, schema_path: str) -> None:
    """Validate the user config, apply schema defaults, and check genome support.

    SalmonTE is disabled with a warning, not an error, for unsupported genomes.

    Renamed 1.x keys listed in DEPRECATED_CONFIG_KEYS are aliased with a
    deprecation warning; every other schema violation (including removed 1.x
    keys) is collected and raised together with a pointer to the migration guide.

    Args:
        config: the Snakemake config dict, modified in place.
        schema_path: path to workflow/schemas/config.schema.yaml.

    Raises:
        WorkflowError: if the config does not satisfy the schema, or enabled
            analyses are unsupported for the configured genome.
    """
    for old, new in DEPRECATED_CONFIG_KEYS.items():
        if old not in config:
            continue
        if new in config:
            raise WorkflowError(
                f"Config has both '{old}' (1.x name) and '{new}'. Remove '{old}'. "
                f"See the migration guide: {MIGRATION_GUIDE}"
            )
        logger.warning(
            f"CONFIG DEPRECATION: '{old}' has been renamed to '{new}'. Please update "
            f"your config file; '{old}' will be removed in a future release. "
            f"See the migration guide: {MIGRATION_GUIDE}"
        )
        config[new] = config.pop(old)

    with open(schema_path, "r") as f:
        schema = yaml.safe_load(f)
    validator_cls = validators.validator_for(schema)
    errors = sorted(
        validator_cls(schema).iter_errors(config),
        key=lambda e: [str(p) for p in e.absolute_path],
    )
    if errors:
        messages = [_describe_config_error(e, config) for e in errors]
        raise WorkflowError(
            "Invalid config file:\n  "
            + "\n  ".join(messages)
            + f"\nSee the 2.0 migration guide: {MIGRATION_GUIDE}"
        )

    apply_schema_defaults(config, schema_path)
    disable_unsupported_salmonTE(config)
    check_genome_support(config)


def giga_to_byte(g):
    g = g - 2  # leave some memory to other processes
    return g * (1024**3)


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
        basename = str(basename).replace(".gz", "")
        basename = Path(basename)

    return basename


def get_sample_name_from_filepath(pattern, test_string):
    m = re.match(pattern, str(test_string))
    if m:
        gd = m.groupdict()
        return gd["sample"]
    else:
        raise ValueError(
            "File path does not match expected pattern: {}".format(test_string)
        )


def get_sample_sheet_path(wildcards):
    sample_sheet_path = next(
        (
            lib["sample_sheet"]
            for lib in config["comparisons"]
            if lib["name"] == wildcards.serie
        ),
        None,
    )
    sample_sheet_path = Path(sample_sheet_path)
    if not sample_sheet_path.is_absolute():
        sample_sheet_path = sample_sheet_path.resolve()
    return sample_sheet_path


def get_samples(wildcards) -> list[str]:
    # We used to extract this from the raw read file patterns via `filepath_pattern`.
    # Now we rely directly on the 'name' column in the sample sheet to allow flexible filenames.
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)
    return sample_sheet["name"].tolist()


def get_samples_names(wildcards) -> list[str]:
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)
    return sample_sheet["name"].tolist()


def get_bw(wildcards):
    """Builds bigwig paths for rule all"""
    o = []
    for lib in config["comparisons"]:
        sample_sheet = pd.read_csv(lib["sample_sheet"])
        samples = sample_sheet["name"].tolist()

        o += expand(
            star_folder.joinpath("{serie}", "{sample}.bw"),
            serie=lib["name"],
            sample=list(samples),
        )
    return o


def get_star_input(wildcards):
    """Builds input paths for STAR alignment testing if a library is single-end or paired-end"""

    samples_names = get_samples_names(wildcards)
    filenames_no_mate = get_samples(wildcards)
    idx = samples_names.index(wildcards.sample)
    s = filenames_no_mate[idx]

    if wildcards.serie in library_names_single:
        ret = trim_reads_folder.joinpath(wildcards.serie, "{0}.fastq.gz".format(s))
    else:
        ret = [
            trim_reads_folder.joinpath(wildcards.serie, "{0}_1.fastq.gz".format(s)),
            trim_reads_folder.joinpath(wildcards.serie, "{0}_2.fastq.gz".format(s)),
        ]
    return ret


def get_params(wildcards, key, default=""):
    """Returns the value of a specific key for the current serie"""
    params = default
    for lib in config["comparisons"]:
        if lib["name"] == wildcards.serie:
            params = lib.get(key, default)
    return params


def get_resolved_param(wildcards, key, nested_key=None, default=None):
    """
    Resolves a parameter by checking:
    1. Comparison-specific override in config["comparisons"]
    2. Global default in config["defaults"]
    3. Provided default fallback
    """
    # Determine the library name (serie)
    serie = getattr(wildcards, "serie", None)

    # Fallback to hash-based lookup if serie is missing
    if serie is None:
        starte_hash = getattr(wildcards, "starte_hash", None)
        if starte_hash and starte_hash in HASH_TO_PARAMS["starTE"]:
            serie = HASH_TO_PARAMS["starTE"][starte_hash].get("serie")

        if serie is None:
            align_hash = getattr(wildcards, "align_hash", None)
            if align_hash and align_hash in HASH_TO_PARAMS["alignment"]:
                serie = HASH_TO_PARAMS["alignment"][align_hash].get("serie")

    # 1. Library override
    if serie:
        for lib in config.get("comparisons", []):
            if lib["name"] == serie:
                val = lib.get(key)
                if nested_key and isinstance(val, dict):
                    val = val.get(nested_key)
                if val is not None:
                    return val

    # 2. Global defaults
    defaults = config.get("defaults", {})
    val = defaults.get(key)
    if nested_key and isinstance(val, dict):
        val = val.get(nested_key)

    if val is not None:
        return val

    # 3. Provided default fallback (last resort)
    if default is not None:
        return default

    # 4. Failure
    param_name = f"{key}.{nested_key}" if nested_key else key
    raise ValueError(
        f"Parameter '{param_name}' could not be resolved for comparison '{serie}'. "
        "Please define it in your config under 'defaults' or within the 'comparisons' entry."
    )


def get_input_size_gb(input_files):
    if isinstance(input_files, str):
        input_files = [input_files]
    total_bytes = 0
    for f in input_files:
        f_path = Path(f)
        if f_path.exists():
            total_bytes += f_path.stat().st_size
    return total_bytes / (1024**3)


def get_star_mem_mb(wildcards, input):
    # Base memory for index (e.g., 32GB)
    base_mem = 32000
    # Extra memory based on input size (1GB of compressed FASTQ ~ 4GB RAM for sorting)
    input_size_gb = get_input_size_gb(input)
    sort_mem = max(5000, input_size_gb * 4000)
    return int(base_mem + sort_mem)


def get_star_runtime(wildcards, input, attempt):
    input_size_gb = get_input_size_gb(input)
    # 1 hour base + 1 hour per GB of compressed FASTQ
    base_runtime = 60
    scaling = input_size_gb * 60
    return int((base_runtime + scaling) * attempt)


def get_fastqc_runtime(wildcards, input, attempt):
    input_size_gb = get_input_size_gb(input)
    # 20 min base + 5 min per GB of input (compressed FASTQ or BAM)
    base_runtime = 20
    scaling = input_size_gb * 5
    return int((base_runtime + scaling) * attempt)


def get_multiqc_runtime(wildcards, input, attempt):
    input_size_gb = get_input_size_gb(input)
    # 20 min base + 10 min per GB of input (FastQC zips + stats files)
    base_runtime = 20
    scaling = input_size_gb * 10
    return int((base_runtime + scaling) * attempt)


def get_featurecounts_mem_mb(wildcards, input):
    input_size_gb = get_input_size_gb(input)
    # Base 4GB + 2GB per GB of input BAM (approximate)
    return int(max(4000, input_size_gb * 2000))


# genome.label prefix -> SalmonTE --reference_genome. SalmonTE ships TE indexes
# for these species only; any other label disables SalmonTE with a warning.
# Keep in sync with the supported list in docs/configuration/index.md.
SALMONTE_GENOMES = {"mm": "mm", "hg": "hs", "dr": "dr", "dm": "dm"}

GTRNADB_URLS = {
    "mm10": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz",
    "mm39": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz",
}


def get_gtrnadb_url(wildcards):
    label = config["genome"]["label"]
    if label in GTRNADB_URLS:
        return GTRNADB_URLS[label]
    raise ValueError(
        f"Unsupported genome label for GtRNAdb: {label}. "
        f"Supported labels are: {', '.join(GTRNADB_URLS.keys())}"
    )


def get_sample_sheet(wildcards):
    """Returns path to sample sheet for current serie"""
    return get_params(wildcards, "sample_sheet")


def get_raw_fasta(wildcards):
    if "fasta_path" in config["genome"]:
        return config["genome"]["fasta_path"]
    else:
        return str(
            references_folder.joinpath(config["genome"]["label"], "raw_fasta.fa")
        )


def get_raw_gtf(wildcards):
    if "gtf_path" in config["genome"]:
        return config["genome"]["gtf_path"]
    else:
        return str(
            references_folder.joinpath(
                config["genome"]["label"], "raw_annotation.gtf.gz"
            )
        )


def get_mimseq_counts(wildcards):
    samples = get_samples_names(wildcards)
    return expand(
        trna_coverage_folder.joinpath(
            "mimseq_{{serie}}", "{sample}", "{sample}_cluster_counts.txt"
        ),
        sample=samples,
    )


def get_trimmomatic_adaptive_input(wildcards):
    """Returns the adaptive params JSON if adaptive trimming is enabled."""
    if wildcards.trim_hash not in HASH_TO_PARAMS["trim"]:
        return []
    params = HASH_TO_PARAMS["trim"][wildcards.trim_hash]
    t_params = params.get("params", {})
    if isinstance(t_params, dict) and t_params.get("adaptive"):
        # We use a path that includes the sample name to avoid collisions if multiple samples share a hash
        # but have slightly different adaptive needs (though they shouldn't if the hash is correct).
        # Actually, samples with the same hash should have the same adaptive params if they are the same data.
        # But FastQC is sample-specific. So we use {sample} in the path.
        return trim_reads_folder.joinpath(
            "_shared", wildcards.trim_hash, f"{wildcards.sample}.adaptive_params.json"
        )
    return []


def get_trimmomatic_command(wildcards):
    """Dynamic derivation of the trimmomatic command string."""
    params = HASH_TO_PARAMS["trim"][wildcards.trim_hash]
    t_params = params.get("params", {})

    # 1. Adaptive mode
    if isinstance(t_params, dict) and t_params.get("adaptive"):
        adaptive_json = trim_reads_folder.joinpath(
            "_shared", wildcards.trim_hash, f"{wildcards.sample}.adaptive_params.json"
        )
        with open(adaptive_json, "r") as f:
            data = json.load(f)
            return data["trimmomatic_params"]

    # 2. Fixed string mode (explicit)
    if isinstance(t_params, str) and t_params != "default":
        return t_params

    # 3. Default fallback
    is_paired = params.get("paired", False)
    suffix = "PE" if is_paired else "SE"
    return f"ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-{suffix}.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"


def parse_filepath(filepath: PosixPath):
    m = re.match(filepath_pattern, str(filepath))
    if m:
        gd = m.groupdict()
        return gd
    raise ValueError("Invalid path: {}".format(filepath))


# Global registries
SAMPLE_HASHES = {}
HASH_TO_PARAMS = {"trim": {}, "alignment": {}, "starTE": {}, "markdup": {}}


def get_sample_hash(serie, sample, step="alignment"):
    """Returns the pre-calculated hash for a given sample and step."""
    return SAMPLE_HASHES.get((serie, sample), {}).get(step, "unknown")


def calculate_content_hash(data_dict):
    """Generates a stable hash from a dictionary of parameters/files."""
    encoded = json.dumps(data_dict, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()[:16]


def populate_sample_registry():
    """
    Scans all sequencing libraries and populates SAMPLE_HASHES and HASH_TO_PARAMS
    for deduplication of trimming and alignment tasks.
    """
    global SAMPLE_HASHES, HASH_TO_PARAMS

    genome_params = {
        "label": config["genome"].get("label"),
        "fasta": config["genome"].get("fasta_path"),
        "gtf": config["genome"].get("gtf_path"),
    }

    for lib in config.get("comparisons", []):
        serie = lib["name"]
        sheet_path = Path(lib["sample_sheet"])
        if not sheet_path.is_absolute():
            sheet_path = sheet_path.resolve()

        if not sheet_path.exists():
            continue

        df = pd.read_csv(sheet_path)
        is_paired = serie in library_names_paired

        raw_trim_params = lib.get("trimmomatic", "default")
        if raw_trim_params == "adaptive":
            trim_params = {"adaptive": True}
        else:
            trim_params = raw_trim_params

        star_params = lib.get("star", "default")

        # New: Resolve starTE and strandedness for hashing
        # We use empty wildcards here since we only care about the 'serie' (lib name)
        w = Wildcards(fromdict={"serie": serie})

        lib_strandedness = get_resolved_param(w, "strandedness", default=0)
        lib_starTE_random = get_resolved_param(w, "starTE_random", default={})
        lib_starTE_multihit = get_resolved_param(w, "starTE_multihit", default={})

        for _, row in df.iterrows():
            sample_name = row["name"]

            # --- Trimming Hash ---
            if is_paired:
                inputs = {
                    "f1": str(resolve_raw_read_path(serie, row.get("filename_1"))),
                    "f2": str(resolve_raw_read_path(serie, row.get("filename_2"))),
                }
            else:
                inputs = {"f": str(resolve_raw_read_path(serie, row.get("filename")))}

            trim_data = {"inputs": inputs, "params": trim_params, "paired": is_paired}
            trim_hash = calculate_content_hash(trim_data)

            # --- Alignment Hash (STAR) ---
            align_data = {
                "trim_hash": trim_hash,
                "star_params": star_params,
                "genome": genome_params,
            }
            align_hash = calculate_content_hash(align_data)

            # --- starTE Hash ---
            starte_data = {
                "trim_hash": trim_hash,
                "genome": genome_params,
                "strandedness": lib_strandedness,
                "starTE_random": lib_starTE_random,
                "starTE_multihit": lib_starTE_multihit,
                "version": "v3_refined_modes",
            }
            starte_hash = calculate_content_hash(starte_data)

            # --- Markdup Hash ---
            markdup_data = {"align_hash": align_hash, "params": "picard_v1_fixed"}
            markdup_hash = calculate_content_hash(markdup_data)

            SAMPLE_HASHES[(serie, sample_name)] = {
                "trim": trim_hash,
                "alignment": align_hash,
                "starTE": starte_hash,
                "markdup": markdup_hash,
            }

            # --- Reverse Mappings ---
            if trim_hash not in HASH_TO_PARAMS["trim"]:
                HASH_TO_PARAMS["trim"][trim_hash] = {
                    "serie": serie,
                    "sample": sample_name,
                    "paired": is_paired,
                }

            if align_hash not in HASH_TO_PARAMS["alignment"]:
                HASH_TO_PARAMS["alignment"][align_hash] = {
                    "serie": serie,
                    "sample": sample_name,
                    "trim_hash": trim_hash,
                }

            if starte_hash not in HASH_TO_PARAMS["starTE"]:
                HASH_TO_PARAMS["starTE"][starte_hash] = {
                    "serie": serie,
                    "sample": sample_name,
                    "paired": is_paired,
                }

            if markdup_hash not in HASH_TO_PARAMS["markdup"]:
                HASH_TO_PARAMS["markdup"][markdup_hash] = {
                    "serie": serie,
                    "sample": sample_name,
                    "align_hash": align_hash,
                }


def get_shared_trim_path(trim_hash, sample, suffix="", ext=".fastq.gz"):
    """Returns the path to a shared trimmed fastq (or an associated stats/summary file)."""
    return trim_reads_folder.joinpath("_shared", trim_hash, f"{sample}{suffix}{ext}")


def get_shared_star_path(align_hash, sample, suffix=""):
    """Returns the path to a shared alignment result."""
    return star_folder.joinpath("_shared", align_hash, f"{sample}{suffix}")


def get_shared_starTE_path(starte_hash, sample, mode="random", suffix=""):
    """Returns the path to a shared starTE result."""
    return starTE_folder.joinpath("_shared", starte_hash, mode, f"{sample}{suffix}")


def get_shared_markdup_path(markdup_hash, sample, suffix=""):
    """Returns the path to a shared markdup result."""
    return markdup_folder.joinpath("_shared", markdup_hash, f"{sample}{suffix}")


def resolve_raw_read_path(serie: str, file_val: str) -> Path:
    """
    Resolves a raw read filename string from the sample sheet into an actual Path.
    - If absolute, returns it directly (testing common extensions if the exact path is missing).
    - If relative, resolves relative to the current working directory.
    - Tries common extensions if they are missing.
    """
    if pd.isna(file_val) or not file_val:
        return None

    p = Path(file_val)

    # Possible extensions if the user didn't provide one
    extensions = ["", ".fastq.gz", ".fq.gz", ".fastq", ".fq"]

    # 1. Check if it's an absolute path
    if p.is_absolute():
        for ext in extensions:
            test_p = p.with_name(p.name + ext)
            if test_p.exists():
                return test_p
        # If absolute path doesn't exist, return as-is so validation fails informatively
        return p

    # 2. Check relative to the current working directory
    for ext in extensions:
        test_p = p.with_name(p.name + ext)
        if test_p.exists():
            return test_p.resolve()

    # If we can't find it, return the CWD-resolved path so snakemake/validation fails informatively
    return p.resolve()


def get_fastq_paired(wildcards):
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)

    # Find the row for this sample
    row = sample_sheet[sample_sheet["name"] == wildcards.sample]
    if row.empty:
        raise ValueError(
            f"Sample {wildcards.sample} not found in sample sheet for serie {wildcards.serie}"
        )

    row = row.iloc[0]

    m1_path = resolve_raw_read_path(wildcards.serie, row.get("filename_1"))
    m2_path = resolve_raw_read_path(wildcards.serie, row.get("filename_2"))

    if not m1_path or not m2_path:
        raise ValueError(
            f"Missing 'filename_1' or 'filename_2' in sample sheet for paired-end sample {wildcards.sample}"
        )

    return {"m1": str(m1_path), "m2": str(m2_path)}


def get_fastq(wildcards):
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)

    # Find the row for this sample
    row = sample_sheet[sample_sheet["name"] == wildcards.sample]
    if row.empty:
        raise ValueError(
            f"Sample {wildcards.sample} not found in sample sheet for serie {wildcards.serie}"
        )

    row = row.iloc[0]

    file_path = resolve_raw_read_path(wildcards.serie, row.get("filename"))

    if not file_path:
        raise ValueError(
            f"Missing 'filename' in sample sheet for single-end sample {wildcards.sample}"
        )

    return str(file_path)


def mkdir(p: Path, verbose=False):
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
        if verbose:
            print("Created {}".format(p))


def get_trna_coverage(wildcards):
    sample_sheet_path = get_sample_sheet_path(wildcards)
    return {
        "bed": expand(
            trna_coverage_folder.joinpath(wildcards.serie, "{sample}.bed"),
            sample=get_samples_names(wildcards),
        ),
        "sample_sheet": sample_sheet_path,
    }


def get_deseq2_test(wildcards):
    return get_resolved_param(wildcards, "deseq2", "test")


def get_deseq2_variable(wildcards):
    return get_resolved_param(wildcards, "deseq2", "variable")


def get_deseq2_reference_level(wildcards):
    return get_resolved_param(wildcards, "deseq2", "reference_level")


def get_salmonTE_quant_input(wildcards):
    sample_sheet_path = get_sample_sheet_path(wildcards)
    sample_sheet = pd.read_csv(sample_sheet_path)

    ret = []
    is_paired = wildcards.serie in library_names_paired

    for _, row in sample_sheet.iterrows():
        if is_paired:
            m1_path = resolve_raw_read_path(wildcards.serie, row.get("filename_1"))
            m2_path = resolve_raw_read_path(wildcards.serie, row.get("filename_2"))
            if m1_path:
                ret.append(m1_path)
            if m2_path:
                ret.append(m2_path)
        else:
            file_path = resolve_raw_read_path(wildcards.serie, row.get("filename"))
            if file_path:
                ret.append(file_path)

    return ret


def set_salmonTE_genome():
    genome = get_salmonTE_genome(config["genome"]["label"])
    if genome is None:
        # unreachable: disable_unsupported_salmonTE() turns SalmonTE off first
        raise ValueError(
            f"Unsupported genome label for salmonTE: {config['genome']['label']}."
        )
    return genome


def get_multiqc_inputs(wildcards):
    """Returns all QC files across every pipeline stage for the consolidated multiqc rule."""
    samples = get_samples_names(wildcards)
    is_paired = wildcards.serie in library_names_paired

    # --- Raw FastQC ---
    if is_paired:
        raw_fastqc = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_1_fastqc.zip"),
            sample=samples,
        ) + expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_2_fastqc.zip"),
            sample=samples,
        )
    else:
        raw_fastqc = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_fastqc.zip"),
            sample=samples,
        )

    # --- Trimmed FastQC + Trimmomatic stats ---
    if is_paired:
        trim_fastqc = expand(
            fastqc_trim_folder.joinpath(
                wildcards.serie, "{sample}_trimmed_1_fastqc.zip"
            ),
            sample=samples,
        ) + expand(
            fastqc_trim_folder.joinpath(
                wildcards.serie, "{sample}_trimmed_2_fastqc.zip"
            ),
            sample=samples,
        )
    else:
        trim_fastqc = expand(
            fastqc_trim_folder.joinpath(wildcards.serie, "{sample}_trimmed_fastqc.zip"),
            sample=samples,
        )

    trimmomatic_stats = expand(
        trim_reads_folder.joinpath(wildcards.serie, "{sample}.stats.txt"),
        sample=samples,
    )

    # --- STAR logs + gene counts + post-alignment FastQC ---
    star_logs = expand(
        star_folder.joinpath(wildcards.serie, "{sample}.Log.final.out"),
        sample=samples,
    )
    star_gene_counts = expand(
        star_folder.joinpath(wildcards.serie, "{sample}.ReadsPerGene.out.tab"),
        sample=samples,
    )
    star_fastqc = expand(
        fastqc_star_folder.joinpath(
            wildcards.serie, "{sample}.Aligned.sortedByCoord.out_fastqc.zip"
        ),
        sample=samples,
    )

    # --- Markdup stats + FastQC ---
    markdup_stats = expand(
        markdup_folder.joinpath(wildcards.serie, "{sample}.markdup.stats.txt"),
        sample=samples,
    )
    markdup_fastqc = expand(
        fastqc_markdup_folder.joinpath(wildcards.serie, "{sample}.markdup_fastqc.zip"),
        sample=samples,
    )

    return {
        "raw_fastqc": raw_fastqc,
        "trim_fastqc": trim_fastqc,
        "trimmomatic_stats": trimmomatic_stats,
        "star_logs": star_logs,
        "star_gene_counts": star_gene_counts,
        "star_fastqc": star_fastqc,
        "markdup_stats": markdup_stats,
        "markdup_fastqc": markdup_fastqc,
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def get_trimmed_fastq(wildcards):
    if wildcards.serie in library_names_single:
        return trim_reads_folder.joinpath(
            wildcards.serie, "{0}.fastq.gz".format(wildcards.sample)
        )
    elif wildcards.serie in library_names_paired:
        return [
            trim_reads_folder.joinpath(
                wildcards.serie, "{0}_1.fastq.gz".format(wildcards.sample)
            ),
            trim_reads_folder.joinpath(
                wildcards.serie, "{0}_2.fastq.gz".format(wildcards.sample)
            ),
        ]
    else:
        raise ValueError(
            f"Could not determine protocol type for comparison '{wildcards.serie}': "
            "set 'protocol: pe' or 'protocol: se' on its entry under 'comparisons' "
            "in the config file."
        )


def get_fastqc_raw(wildcards):
    """Returns raw fastqc zip paths; used as input tracker in fastqc_raw rules."""
    s = get_samples(wildcards)
    if wildcards.serie in library_names_single:
        fastqcs = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_fastqc.zip"),
            sample=s,
        )
    else:
        fastqcs = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_1_fastqc.zip"),
            sample=s,
        ) + expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_raw_2_fastqc.zip"),
            sample=s,
        )
    return {"fastqc": fastqcs, "sample_sheet": get_sample_sheet_path(wildcards)}


def build_rule_all_inputs(wildcards):
    ret = []

    # MultiQC reports
    ret += expand(
        multiqc_folder.joinpath("{serie}", "multiqc_report.html"),
        serie=library_names_single + library_names_paired,
    )

    # DESeq2 flags
    ret += expand(
        rdata_folder.joinpath("deseq2/{serie}/dds.rds"),
        serie=library_names_single + library_names_paired,
    )
    ret += expand(
        analysis_folder.joinpath("datavzrd", "{serie}", "datavzrd"),
        serie=library_names_single + library_names_paired,
    )
    # Bigwig files
    ret += get_bw(wildcards)

    # TE
    if not config["disable_TE_analysis"]:
        # FeatureCounts tables from STAR-TE
        ret += expand(
            starTE_folder.joinpath("{se_serie}/featureCount/{method}.txt"),
            se_serie=library_names_single,
            method=["multihit", "random"],
        )
        ret += expand(
            starTE_folder.joinpath("{pe_serie}/featureCount/{method}.txt"),
            pe_serie=library_names_paired,
            method=["multihit", "random"],
        )
        ret += expand(
            starTE_folder.joinpath("{serie}", "random", "datavzrd"),
            serie=library_names_single + library_names_paired,
        )

    # tRNA coverage files
    if not config["disable_tRNA_analysis"]:
        method = config.get("tRNA_quantification", {}).get("method", "standard")
        m_str = "mimseq" if method == "mim-tRNA-seq" else "standard"

        ret += expand(
            trna_coverage_folder.joinpath("{serie}", f"tRNA_lfc_{m_str}.csv"),
            serie=library_names_paired + library_names_single,
        )
        ret += expand(
            trna_coverage_folder.joinpath("{serie}", f"datavzrd_{m_str}"),
            serie=library_names_paired + library_names_single,
        )

    # SalmonTE quantification
    if not config["disable_salmonTE_analysis"]:
        # SalmonTE quantification result
        ret += expand(
            data_folder.joinpath("salmonTE/de_analysis/{se_serie}"),
            se_serie=library_names_single,
        )
        ret += expand(
            data_folder.joinpath("salmonTE/de_analysis/{pe_serie}"),
            pe_serie=library_names_paired,
        )

    return ret


def get_deseq2_inputs(wildcards):
    return {
        "star_counts": expand(
            star_folder.joinpath(wildcards.serie, "{sample}.ReadsPerGene.out.tab"),
            sample=get_samples_names(wildcards),
        ),
        "annotation_file": gtf_path,
        "sample_sheet": get_sample_sheet_path(wildcards),
    }


def validate_sample_sheets():
    """Validates that all specified raw reads files in sample sheets exist."""
    print("Validating sample sheets against filesystem...")
    missing_files = []

    for lib in config.get("comparisons", []):
        serie = lib["name"]
        sheet_path = Path(lib["sample_sheet"])
        # If relative, it's relative to CWD
        if not sheet_path.is_absolute():
            sheet_path = sheet_path.resolve()

        if not sheet_path.exists():
            missing_files.append(f"Sample sheet missing: {sheet_path}")
            continue

        df = pd.read_csv(sheet_path)

        is_paired = serie in library_names_paired

        required_columns = (
            ["name", "filename_1", "filename_2"] if is_paired else ["name", "filename"]
        )
        missing_columns = [c for c in required_columns if c not in df.columns]
        if missing_columns:
            missing_files.append(
                f"Sample sheet {sheet_path} lacks column(s) {', '.join(missing_columns)} "
                f"required by 'protocol: {lib.get('protocol')}' of comparison '{serie}' "
                "in the config file."
            )
            continue

        for _, row in df.iterrows():
            sample_name = row["name"]
            if is_paired:
                m1 = resolve_raw_read_path(serie, row.get("filename_1"))
                m2 = resolve_raw_read_path(serie, row.get("filename_2"))
                if m1 is None or not Path(m1).exists():
                    missing_files.append(
                        f"Serie '{serie}' / Sample '{sample_name}': Could not resolve filename_1 ('{row.get('filename_1')}')"
                    )
                if m2 is None or not Path(m2).exists():
                    missing_files.append(
                        f"Serie '{serie}' / Sample '{sample_name}': Could not resolve filename_2 ('{row.get('filename_2')}')"
                    )
            else:
                f = resolve_raw_read_path(serie, row.get("filename"))
                if f is None or not Path(f).exists():
                    missing_files.append(
                        f"Serie '{serie}' / Sample '{sample_name}': Could not resolve filename ('{row.get('filename')}')"
                    )

    if missing_files:
        error_msg = "\n" + "\n".join(missing_files)
        raise WorkflowError(f"Sample sheet validation failed:{error_msg}")
# We don't call it here anymore because library_names_paired isn't populated until Snakefile line 110.
# It is called from the Snakefile.
