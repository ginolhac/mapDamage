import unittest
import subprocess


def read_nocomment(dnacomp):
    with open(dnacomp, 'r') as f:
        a = f.readlines()
        return [x for x in a if not x.startswith('#')]

class testCases(unittest.TestCase):

    def test_no_stats(self):
        """regular mapDamage"""
        subprocess.run(["mapDamage",  "-i",  "tests/test.bam", "-r", "tests/fake1.fasta", "-d", "tests/results", "--no-stats"], check = True)
        self.assertTrue(read_nocomment("tests/dnacomp.txt") == read_nocomment("tests/results/dnacomp.txt"))

if  __name__=='__main__':
    unittest.main()
