import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_download_repeatmasker_annotation_file():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(
            ".tests/unit/download_repeatmasker_annotation_file/data"
        )
        expected_path = PurePosixPath(
            ".tests/unit/download_repeatmasker_annotation_file/expected"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print("references/rmsk/mm10.gtf references/rmsk/mm10.bed", file=sys.stderr)

        # Run the test job
        cmd = utils.get_cmd(
            workdir, ["references/rmsk/mm10.gtf", "references/rmsk/mm10.bed"], cwd
        )
        sp.check_output(cmd)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check()
