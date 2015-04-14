import os
import unittest

import analysis.bed_pattern_matching as bpm


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bed_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/in.bed")
        cls.fasta_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/fasta.fa")
        cls.pfm_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/Sox2.dat")

        cls.parser = bpm.create_parser()

    def test_full_args(self):
        args = self.parser.parse_args([
            str(self.fasta_file_path),
            str(self.pfm_file_path),
            str(self.bed_file_path),
            "-tf", "Sox2",
            "-th", "0.8",
            "-rc",
            "-c"
        ])


if __name__ == "__main__":
    unittest.main()