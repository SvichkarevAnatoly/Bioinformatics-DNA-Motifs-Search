import unittest
from src.lib import parse_interval_line


class Test(unittest.TestCase):
    def test_parse_interval_line(self):
        test_interval_str = "chr1:5-25"
        actual_interval = parse_interval_line(test_interval_str)

        expected_interval = ["chr1", 5, 25]

        self.assertEqual(expected_interval, actual_interval)


if __name__ == "__main__":
    unittest.main()