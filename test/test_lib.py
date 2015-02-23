from Bio import motifs
import os
import unittest
import src.lib as lib


class Test(unittest.TestCase):
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

    def test_create_matrix_from_pwms(self):
        # TODO: don't use reading file
        pwms_filename = os.path.join(os.path.dirname(__file__), "test_data/pwms_transfac.dat")
        with open(pwms_filename, 'r') as pwms_file:
            pwm_record_list = motifs.parse(pwms_file, "TRANSFAC")

        actual_matrix = lib.create_matrix_from_pwms("motif1", pwm_record_list)
        expected_matrix = [
            [1.0, 2.0, 3.0],
            [1.0, 1.0, 4.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0]
        ]   # N    A    M

        self.assertEqual(expected_matrix, actual_matrix)


if __name__ == "__main__":
    unittest.main()