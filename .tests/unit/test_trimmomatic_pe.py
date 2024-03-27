import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_trimmomatic_pe():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/trimmomatic_pe/data")
        expected_path = PurePosixPath(".tests/unit/trimmomatic_pe/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print(
            "results/trim/GSE130735-subset/SRX5795113_SRR9016959_1.fastq.gz results/trim/GSE130735-subset/SRX5795113_SRR9016959_2.fastq.gz results/trim/GSE130735-subset/SRX5795113_SRR9016959_1.unpaired.fastq.gz results/trim/GSE130735-subset/SRX5795113_SRR9016959_2.unpaired.fastq.gz results/trim/GSE130735-subset/SRX5795113_SRR9016959.summary.txt",
            file=sys.stderr,
        )

        # Run the test job.
        cmd = utils.get_cmd(
            workdir,
            [
                "results/trim/GSE130735-subset/SRX5795113_SRR9016959_1.fastq.gz",
                "results/trim/GSE130735-subset/SRX5795113_SRR9016959_2.fastq.gz",
                "results/trim/GSE130735-subset/SRX5795113_SRR9016959_1.unpaired.fastq.gz",
                "results/trim/GSE130735-subset/SRX5795113_SRR9016959_2.unpaired.fastq.gz",
                "results/trim/GSE130735-subset/SRX5795113_SRR9016959.summary.txt",
            ],
            cwd,
        )
        sp.check_output(cmd)

        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "-f",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--configfile",
                ".tests/integration/config.yaml",
                "--use-conda",
                "--use-apptainer",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check(
            ignore_pattern=r".*\.stats.txt"
        )
