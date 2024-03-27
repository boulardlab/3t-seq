import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_coverage_trna():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/coverage_trna/data")
        expected_path = PurePosixPath(".tests/unit/coverage_trna/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print(
            "results/tRNA_coverage/GSE130735-subset/SRX5795117_SRR9016963.bed",
            file=sys.stderr,
        )

        # Run the test job.
        cmd = utils.get_cmd(
            workdir,
            "results/tRNA_coverage/GSE130735-subset/SRX5795117_SRR9016963.bed",
            cwd,
        )
        sp.check_output(cmd)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from utils.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check()
