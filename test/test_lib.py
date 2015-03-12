import unittest
import cStringIO

from Bio import motifs

import src.lib as lib


class TestLib(unittest.TestCase):
    @classmethod
    def get_pwm_records(cls):
        pwms_str = '\n'.join([
            "VV  January 28, 2015 08:46:03",
            "XX", "//", "ID  motif1",
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
        ])

        tempfile = cStringIO.StringIO()
        tempfile.write(pwms_str)
        tempfile.seek(0)
        pwm_records = motifs.parse(tempfile, "TRANSFAC")
        tempfile.close()
        return pwm_records

    @classmethod
    def setUpClass(cls):
        cls.pwm_records = cls.get_pwm_records()

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
            [0, 0]  # T
        ]
        matching = lib.search_motif(seq, [matrix], 0.7, False)
        matching = matching[0][0]

        expected_position = 1
        self.assertEqual(expected_position, matching[0])
        expected_score = 2.0
        self.assertEqual(expected_score, matching[1])

    def test_search_motif_both_strand(self):
        seq = "ACGTGC"
        matrix = [
            [0, 0],  # A
            [1, 0],  # C
            [0, 1],  # G
            [0, 0]  # T
        ]
        results = lib.search_motif(seq, [matrix], 0.7, True)[0]
        pos1 = results[0]
        pos2 = results[1]

        expected_pos1 = (1, 2.0)
        self.assertEqual(expected_pos1, pos1)
        expected_pos2 = (-5, 2.0)
        self.assertEqual(expected_pos2, pos2)

    def test_search_motif_two_matrices(self):
        seq = "ACGTGC"
        matrix1 = [
            [0, 0, 0],  # A
            [0, 0, 0],  # C
            [1, 0, 1],  # G
            [0, 1, 0]  # T
        ]
        matrix2 = [
            [1, 0, 0],  # A
            [0, 1, 0],  # C
            [0, 0, 1],  # G
            [0, 0, 0]  # T
        ]
        results = lib.search_motif(seq, [matrix1, matrix2], 0.7, False)
        matrix1_pos = results[0][0]
        matrix2_pos = results[1][0]

        expected_matrix1_pos = (2, 3.0)
        self.assertEqual(expected_matrix1_pos, matrix1_pos)
        expected_matrix2_pos = (0, 3.0)
        self.assertEqual(expected_matrix2_pos, matrix2_pos)

    def test_get_pwm_ids(self):
        actual_pwm_ids = lib.get_pwm_ids(self.pwm_records)
        expected_pwm_ids = ["MOTIF1", "MOTIF2"]
        self.assertEquals(expected_pwm_ids, actual_pwm_ids)

    def test_filter_pwms_in_tfs(self):
        tf_names = ["MOTIF1"]
        actual_pwms = lib.filter_pwms_in_tfs(self.pwm_records, tf_names)
        expected_pwms = [self.pwm_records[0]]
        self.assertEquals(expected_pwms, actual_pwms)

    def test_filter_tfs_in_pwms_throw_exception_if_tf_not_in_pwms(self):
        tf_names = ["motif1", "motifNotInPwms"]
        pwm_ids = ["motif1", "motif2"]

        # TODO:
        with self.assertRaises(Exception):
            lib.filter_pwms_in_tfs(tf_names, pwm_ids)

    def test_create_matrices_from_pwms(self):
        tf_names = ["MOTIF1"]
        actual_pwms = lib.create_matrices_from_pwms(self.pwm_records, tf_names)
        expected_matrix1 = [
            [1.0, 2.0, 3.0],
            [1.0, 1.0, 4.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0]
        ]  # N    A    M
        expected_pwms = [expected_matrix1]
        self.assertEqual(expected_pwms, actual_pwms)

    def test_get_join_position_str_positive_positions(self):
        positions = [1, 2, 3]
        seq_length = 10

        actual_position_str = lib.get_join_position_str(positions, seq_length)
        expected_position_str = "1;2;3"

        self.assertEquals(expected_position_str, actual_position_str)

if __name__ == "__main__":
    unittest.main()