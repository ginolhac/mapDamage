import filecmp
import optparse
import subprocess
import unittest

import pysam

import mapdamage


def read_nocomment(dnacomp):
    with open(dnacomp, "r") as f:
        a = f.readlines()
        return [x for x in a if not x.startswith("#")]


def mock_options(filename, rescale_out, folder):
    """Make the options object with nice values for testing"""
    return optparse.Values(
        {
            "filename": filename,
            "rescale_out": rescale_out,
            "verbose": True,
            "folder": folder,
            "rescale_length_5p": 12,  # default values as in --seq-length
            "rescale_length_3p": 12,  # default values as in --seq-length
            "quiet": True,
        }
    )


class testCases(unittest.TestCase):
    def test_no_stats(self):
        """regular mapDamage"""
        subprocess.run(
            [
                "mapDamage",
                "-i",
                "tests/test.bam",
                "-r",
                "tests/fake1.fasta",
                "-d",
                "tests/results",
                "--no-stats",
            ],
            check=True,
        )
        self.assertTrue(
            read_nocomment("tests/dnacomp.txt")
            == read_nocomment("tests/results/dnacomp.txt")
        )


class testRescaling(unittest.TestCase):
    def test_single_end_file(self):
        """Test, rescaling BAM file"""
        #
        # The expected substitution frequencies before and after scaling using the scaled qualities as probalities:
        # CT	0.06226411977920493		0.04163524443356556
        # TC	0.020395286584806528		0.020395286584806528
        # GA	0.04400459948304954		0.03794905109091021
        # AG	0.05355350777642983		0.05355350777642983
        # Quality metrics before and after scaling
        # CT-Q0 	5		5
        # CT-Q10 	5		2
        # CT-Q20 	3		2
        # CT-Q30 	3		2
        # CT-Q40 	0		0
        # GA-Q0 	5		5
        # GA-Q10 	5		4
        # GA-Q20 	1		0
        # GA-Q30 	1		0
        # GA-Q40 	0		0
        options = mock_options(
            "tests/test.bam", "tests/test.rescaled.sam", "tests/probs/"
        )
        ref = pysam.Fastafile("tests/fake1.fasta")
        mapdamage.rescale.rescale_qual(ref, options, debug=True)
        self.assertTrue(
            filecmp.cmp("tests/test.rescaled.sam", "tests/test.rescaled.correct.sam")
        )


if __name__ == "__main__":
    unittest.main()
