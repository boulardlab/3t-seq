import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_deseq2_tRNA():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/deseq2_tRNA/data")
        expected_path = PurePosixPath(".tests/unit/deseq2_tRNA/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print(
            "results/tRNA_coverage/GSE130735-subset/tRNA_dds.rds results/tRNA_coverage/GSE130735-subset/tRNA_lfc.txt",
            file=sys.stderr,
        )

        # Run the test job.
        cmd = utils.get_cmd(
            workdir,
            [
                "results/tRNA_coverage/GSE130735-subset/tRNA_dds.rds",
                "results/tRNA_coverage/GSE130735-subset/tRNA_lfc.txt",
            ],
            cwd,
        )
        sp.check_output(cmd)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check()
