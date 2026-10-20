"""Build Trimmomatic step strings from base steps and user `extra_params`.

Used by workflow/rules/common.smk (non-adaptive `extra_params`) and by
workflow/scripts/derive_trim_params.py (adaptive mode), so both forms of the
`trimmomatic` setting merge user steps the same way (issue #48).
"""

ADAPTERS_DIR = "$CONDA_PREFIX/share/trimmomatic/adapters"


def step_type(step: str) -> str:
    """Return the step name, e.g. "MINLEN" for "MINLEN:36"."""
    return step.split(":", 1)[0].upper()


def default_steps(is_paired: bool) -> list[str]:
    """Built-in default Trimmomatic steps."""
    suffix = "PE" if is_paired else "SE"
    return [
        f"ILLUMINACLIP:{ADAPTERS_DIR}/TruSeq3-{suffix}.fa:2:30:10",
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        "MINLEN:36",
    ]


def merge_steps(base_steps: list[str], extra_params: str = "") -> str:
    """Merge user steps into base steps and return the Trimmomatic step string.

    Rules:
    - a user step replaces every base step of the same type;
    - order: ILLUMINACLIP, remaining base steps, user steps, MINLEN;
    - if CROP:N is present and MINLEN comes from the base steps, MINLEN is
      capped at N, so cropping cannot make every read fail MINLEN.

    Args:
        base_steps: derived or default steps, in their intended order.
        extra_params: whitespace-separated user steps, may be empty.

    Returns:
        Space-separated Trimmomatic steps.

    Raises:
        ValueError: if a CROP or MINLEN value is not an integer.
    """
    user_steps = (extra_params or "").split()
    user_types = {step_type(s) for s in user_steps}
    base_kept = [s for s in base_steps if step_type(s) not in user_types]

    def pick(kind: str) -> list[str]:
        return [s for s in user_steps + base_kept if step_type(s) == kind]

    clip = pick("ILLUMINACLIP")
    minlen = pick("MINLEN")
    middle = [s for s in base_kept if step_type(s) not in ("ILLUMINACLIP", "MINLEN")]
    middle += [s for s in user_steps if step_type(s) not in ("ILLUMINACLIP", "MINLEN")]

    crop = [s for s in user_steps + base_kept if step_type(s) == "CROP"]
    if crop and minlen and "MINLEN" not in user_types:
        crop_len = int(crop[-1].split(":")[1])
        minlen = [f"MINLEN:{min(int(m.split(':')[1]), crop_len)}" for m in minlen]

    return " ".join(clip + middle + minlen)
