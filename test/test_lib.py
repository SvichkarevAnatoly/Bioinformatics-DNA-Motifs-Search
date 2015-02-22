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

    def test_extend_interval(self):
        test_interval_param = ["chr1", 55, 75]
        actual_interval = lib.interval_extend(test_interval_param, 100)

        expected_interval = ["chr1", 15, 115]

        self.assertEqual(expected_interval, actual_interval)


if __name__ == "__main__":
    unittest.main()