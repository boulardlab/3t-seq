import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import utils


def test_star():

    cwd = Path().resolve()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/star/data")
        expected_path = PurePosixPath(".tests/unit/star/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        utils.prepenv(workdir)

        # dbg
        print(
            "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Aligned.sortedByCoord.out.bam results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Aligned.toTranscriptome.out.bam results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.ReadsPerGene.out.tab results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.SJ.out.tab results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.Unique.str1.out.wig results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.Unique.str2.out.wig results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.UniqueMultiple.str1.out.wig results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.UniqueMultiple.str2.out.wig",
            file=sys.stderr,
        )

        # Run the test job.
        cmd = utils.get_cmd(
            workdir,
            [
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Aligned.sortedByCoord.out.bam",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Aligned.toTranscriptome.out.bam",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.ReadsPerGene.out.tab",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.SJ.out.tab",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.Unique.str1.out.wig",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.Unique.str2.out.wig",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.UniqueMultiple.str1.out.wig",
                "results/alignments/star/GSE130735-subset/SRX5795117_SRR9016963.Signal.UniqueMultiple.str2.out.wig",
            ],
            cwd,
        )
        sp.check_output(cmd)

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        utils.OutputChecker(data_path, expected_path, workdir).check()
