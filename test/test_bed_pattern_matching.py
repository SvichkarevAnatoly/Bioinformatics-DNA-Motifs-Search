import os
import unittest
import sys
import cStringIO

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
            "-c", "100"
        ])
        self.assertEqual(sys.stdout, args.output)
        self.assertEqual(["SOX2"], args.tf)
        self.assertEqual(0.8, args.threshold)
        self.assertTrue(args.reverse_complement)
        self.assertEquals(100, args.constriction)

    def test_read_bed_file(self):
        args = self.parser.parse_args([
            str(self.fasta_file_path),
            str(self.pfm_file_path),
            str(self.bed_file_path)
        ])

        bed_peaks = args.bed
        expected_peak1 = ("chr1", "4736010", "4736158", "Z4_Sox2_peak_1")
        self.assertSequenceEqual(expected_peak1, bed_peaks.next())

        expected_peak2 = ("chr1", "5223047", "5223196", "Z4_Sox2_peak_2")
        self.assertSequenceEqual(expected_peak2, bed_peaks.next())

    def test_example_from_data_samples(self):
        args = self.parser.parse_args([
            str(self.fasta_file_path),
            str(self.pfm_file_path),
            str(self.bed_file_path),
            "-tf", "Sox2",
            "-th", "0.7",
            "-rc",
            "-c", "100"
        ])
        args.output = cStringIO.StringIO()

        result = bpm.process(args)
        bpm.save(result, args)

        args.output.seek(0)
        actual_file_contents = args.output.read()
        expected_contents = "\n".join([
            "[Z4_Sox2_peak_1]\t56;85;119;4(-);111(-);128(-)\tCCTGTGCT",
            "[Z4_Sox2_peak_2]\t14;27;41;6(-);94(-);107(-);126(-)\tCCTTTGTT"
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

if __name__ == "__main__":
    unittest.main()