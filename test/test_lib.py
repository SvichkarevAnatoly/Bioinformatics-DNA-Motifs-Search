import unittest
import cStringIO

from Bio import motifs

import src.lib as lib


class TestLib(unittest.TestCase):
    # TODO: test reading pws
    def test_parse_interval_line(self):
        test_interval_str = "chr1:5-25"
        actual_interval = lib.parse_interval_line(test_interval_str)

        expected_interval = ["chr1", 5, 25]

        self.assertEqual(expected_interval, actual_interval)

    def test_interval_param_list_to_str(self):
        test_interval_param = ["chr1", 5, 25]
        actual_interval = lib.interval_param_list_to_str(test_interval_param)

        expected_interval = "chr1:5-25"

        self.assertEqual(expected_interval, actual_interval)

    def test_interval_param_to_str(self):
        chr_name = "chr1"
        start = 5
        end = 25
        actual_interval = lib.interval_param_to_str(chr_name, start, end)

        expected_interval = "chr1:5-25"

        self.assertEqual(expected_interval, actual_interval)

    def test_interval_extend(self):
        new_length = 100
        expected_interval = ["chr1", 15, 115]

        test_interval_param1 = ["chr1", 55, 75]
        actual_interval1 = lib.interval_extend(test_interval_param1, new_length)
        self.assertEqual(expected_interval, actual_interval1)

        test_interval_param2 = ["chr1", 56, 75]
        actual_interval2 = lib.interval_extend(test_interval_param2, new_length)
        self.assertEqual(expected_interval, actual_interval2)

        test_interval_param3 = ["chr1", 55, 76]
        actual_interval3 = lib.interval_extend(test_interval_param3, new_length)
        self.assertEqual(expected_interval, actual_interval3)

    def test_interval_length(self):
        test_interval_param = ["chr1", 55, 75]
        actual_length = lib.interval_length(test_interval_param)

        expected_length = 20

        self.assertEqual(expected_length, actual_length)

    def test_create_output_file_name(self):
        input_file_name = "/home/input.txt"
        actual_file_name = lib.create_output_file_name(input_file_name)

        expected_file_name = "/home/input_out.txt"

        self.assertEqual(expected_file_name, actual_file_name)

    def test_search_motif(self):
        seq = "ACGT"
        matrix = [
            [0, 0],  # A
            [1, 0],  # C
            [0, 1],  # G
            [0, 0]   # T
        ]
        matching = lib.search_motif(seq, [matrix], 0.7, False)
        matching = matching[0][0]

        expected_position = 1
        self.assertEqual(expected_position, matching[0])
        expected_score = 2.0
        self.assertEqual(expected_score, matching[1])

    def test_create_matrices_from_pwms(self):
        tempfile = cStringIO.StringIO()
        tempfile.write('\n'.join([
            "VV  January 28, 2015 08:46:03",
            "XX",
            "//",
            "ID  motif1",
            "P0      A      C      G      T",
            "01      1      1      1      2      N",
            "02      2      1      0      1      A",
            "03      3      4      0      0      M",
            "XX",
            "//",
            "ID  motif2",
            "P0      A      C      G      T",
            "01      1      1      0      0      M",
            "02      2      0      1      0      A",
            "XX",
            "//"
        ]))
        tempfile.seek(0)

        pwm_record_list = motifs.parse(tempfile, "TRANSFAC")

        actual_matrix1 = lib.create_matrices_from_pwms(pwm_record_list, "motif1")
        expected_matrix1 = [
            [1.0, 2.0, 3.0],
            [1.0, 1.0, 4.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0]
        ]  # N    A    M
        self.assertEqual([expected_matrix1], actual_matrix1)

        actual_matrix2 = lib.create_matrices_from_pwms(pwm_record_list, None)
        expected_matrix2 = [
            expected_matrix1,
            [
                [1.0, 2.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [0.0, 0.0]
            ]  # M    A
        ]
        self.assertEqual(expected_matrix2, actual_matrix2)


if __name__ == "__main__":
    unittest.main()