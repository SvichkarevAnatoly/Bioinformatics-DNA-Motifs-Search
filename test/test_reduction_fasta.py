import os
import unittest
import cStringIO

import analysis.reduction_fasta as rf


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bed_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/in.bed")
        cls.raw_fasta_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/raw_fasta.fa")

        cls.parser = rf.create_parser()

    def test_on_raw_fasta_example(self):
        args = self.parser.parse_args([
            str(self.raw_fasta_file_path),
            str(self.bed_file_path),
        ])
        args.output = cStringIO.StringIO()

        rf.save(args)

        args.output.seek(0)
        actual_file_contents = args.output.read()
        expected_contents = "\n".join([
            ">Z4_Sox2_peak_1",
            "GGATAATAAATGACGTGTCCACATGCATTACTTTAGTAAGGTGCAATGCCTTGACGCCTG",
            "TGCTTGTACTAACAGATTTCAACAGCAATTCTTCTTGAATTCCTTGAGTTAATCAATAGC",
            "CTTATTTAAATACTGGAAACTTACTTTT",
            ">Z4_Sox2_peak_2",
            "CTAAGCAACACTGACTTCTGTTTTCCCCCTTTGTTCTGTTTCCTTTGCTCTATACAGCTC",
            "AAAAGAAGCTCCTTTTGAGTCAGGCACAGCAGCAAACAAAAGCAAGTGACAAAGGAGTTG",
            "AAGATAACCAGAGGGCTCCCAGTCTGGGG"
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)


if __name__ == "__main__":
    unittest.main()