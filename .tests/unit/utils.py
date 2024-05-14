from pathlib import Path
import subprocess as sp
import os
import shutil
import yaml
import re
import struct


from common import OutputChecker as BasicOutpuChecker


class OutputChecker(BasicOutpuChecker):
    def __init__(self, data_path, expected_path, workdir):
        super().__init__(data_path, expected_path, workdir)

    def get_input_files(self):
        return set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )

    def get_expected_files(self):
        return set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )

    def check(self, ignore_pattern=None):

        input_files = self.get_input_files()
        expected_files = self.get_expected_files()

        unexpected_files = set()

        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)

                if (
                    str(f).startswith(".snakemake")
                    or str(f).startswith("config.yaml")
                    or str(f).startswith("sample-sheet.csv")
                    or str(f).endswith(".log")
                    or "resources" in str(f.resolve())
                    or (ignore_pattern and re.match(ignore_pattern, str(f)))
                ):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)

        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_bam_files(self, file1_path, file2_path):
        # Open the BAM files
        with open(file1_path, "rb") as file1, open(file2_path, "rb") as file2:
            while True:
                # Read the next 4 bytes to get the length of the record
                len_bytes1 = file1.read(4)
                len_bytes2 = file2.read(4)

                if not len_bytes1 and not len_bytes2:
                    # End of both files
                    return 0
                elif len_bytes1 and len_bytes2:
                    # Unpack the length bytes as an unsigned integer in little-endian byte order
                    record_len1 = struct.unpack("<I", len_bytes1)[0]
                    record_len2 = struct.unpack("<I", len_bytes2)[0]

                    # Read the record data
                    record_data1 = file1.read(record_len1 - 4)
                    record_data2 = file2.read(record_len2 - 4)

                    # Compare the records excluding the header
                    if record_data1 != record_data2:
                        return 1

                else:
                    return 1

    def compare_files(self, generated_file, expected_file):
        if generated_file.suffix == ".bam" and expected_file.suffix == ".bam":
            self.compare_bam_files(generated_file, expected_file)
        else:
            sp.check_output(["cmp", generated_file, expected_file])


def prepenv(workdir):

    configfile_path = ".tests/integration/config.yaml"
    shutil.copy2(configfile_path, workdir)

    with open(configfile_path, "r") as fh:
        configfile = yaml.safe_load(fh)

        for sequencing_library in configfile["sequencing_libraries"]:
            sample_sheet_path = Path(sequencing_library["sample_sheet"])

            new_sample_sheet_path = workdir / sample_sheet_path.parent
            new_sample_sheet_path.mkdir(parents=True, exist_ok=True)
            shutil.copy(sample_sheet_path, new_sample_sheet_path)


def get_cmd(workdir, expected_output, project_root):

    ret = [
        "python",
        "-m",
        "snakemake",
        None,
        "-f",
        "-j1",
        "--target-files-omit-workdir-adjustment",
        "--configfile",
        ".tests/integration/config.yaml",
        "--sdm",
        "conda",
        "apptainer",
        "--apptainer-args",
        "-B {}".format(str(project_root)),
        "--directory",
        workdir,
    ]

    if isinstance(expected_output, list):
        ret.pop(3)
        while expected_output:
            elm = expected_output.pop()
            ret.insert(3, elm)
    else:
        ret[3] = expected_output

    return ret
