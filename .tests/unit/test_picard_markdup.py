import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_picard_markdup():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/picard_markdup/data")
        expected_path = PurePosixPath(".tests/unit/picard_markdup/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print(
            "results/alignments/star_markdup/GSE130735-subset/SRX5795117_SRR9016963.markdup.bam results/alignments/star_markdup/GSE130735-subset/SRX5795117_SRR9016963.markdup.stats.txt",
            file=sys.stderr,
        )

        # Run the test job.
        cmd = utils.get_cmd(
            workdir,
            [
                "results/alignments/star_markdup/GSE130735-subset/SRX5795117_SRR9016963.markdup.bam",
                "results/alignments/star_markdup/GSE130735-subset/SRX5795117_SRR9016963.markdup.stats.txt",
            ],
            cwd,
        )
        sp.check_output(cmd)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check()
