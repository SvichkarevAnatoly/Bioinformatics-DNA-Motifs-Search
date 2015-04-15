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
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.71678019764\t+\t56\t64\t0\tCCTGTGCT",
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.70354024753\t+\t85\t93\t0\tCAATTCTT",
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.71678019764\t+\t119\t127\t0\tCCTTATTT",
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.71678019764\t-\t4\t12\t0\tCATTTATT",
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.70354024753\t-\t111\t119\t0\tCTATTGAT",
            "chr1\t4736010\t4736158\tZ4_Sox2_peak_1\t0.70354024753\t-\t128\t136\t0\tCCAGTATT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.71678019764\t+\t14\t22\t0\tCTTCTGTT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t1.00000000000\t+\t27\t35\t0\tCCTTTGTT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.85839009882\t+\t41\t49\t0\tCCTTTGCT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.70354024753\t-\t6\t14\t0\tTCAGTGTT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.85839009882\t-\t94\t102\t0\tCTTTTGTT",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.96728389139\t-\t107\t115\t0\tCCTTTGTC",
            "chr1\t5223047\t5223196\tZ4_Sox2_peak_2\t0.71678019764\t-\t126\t134\t0\tCCTCTGGT"
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

    def test_seq_search_results_best_matches(self):
        tf_name = "NANOG"
        seq_search_result = bpm.SeqSearchResults("", "", None, [tf_name], [3])
        matches = [(3, 4.0), (1, 6.0), (-6, 3.0)]
        seq_search_result.fill_matches([matches], [6])

        actual_best_match = seq_search_result.best_matches(tf_name)
        expected_best_match = [matches[1]]
        self.assertEqual(expected_best_match, actual_best_match)

    def test_seq_search_results_nearest_to_center_best_matches(self):
        tf_name = "NANOG"
        sequence = "AAAAAAAAAA"
        seq_search_result = bpm.SeqSearchResults("", sequence, None, [tf_name], [3])
        matches = [(3, 4.0), (1, 6.0), (-6, 3.0), (4, 6.0)]
        seq_search_result.fill_matches([matches], [6])

        actual_best_match = seq_search_result.best_matches(tf_name)
        expected_best_match = [matches[3]]
        self.assertEqual(expected_best_match, actual_best_match)

if __name__ == "__main__":
    unittest.main()