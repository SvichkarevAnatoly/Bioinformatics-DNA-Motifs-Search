import unittest
import src.lib as lib


class Test(unittest.TestCase):
    def test_parse_interval_line(self):
        test_interval_str = "chr1:5-25"
        actual_interval = lib.parse_interval_line(test_interval_str)

        expected_interval = ["chr1", 5, 25]

        self.assertEqual(expected_interval, actual_interval)

    def test_interval_param_to_str(self):
        test_interval_param = ["chr1", 5, 25]
        actual_interval = lib.interval_param_to_str(test_interval_param)

        expected_interval = "chr1:5-25"

        self.assertEqual(expected_interval, actual_interval)


if __name__ == "__main__":
    unittest.main()