import os
from strop import strip
import unittest

import src.bed_center_extender as bce


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), 'test_data/test_bed_file.bed')


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = bce.create_parser()

    def setUp(self):
        self.test_data = open(TEST_DATA_FILENAME, 'r')

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_only_bedfile_cl_args(self):
        args = self.parser.parse_args(["test/test_data/test_bed_file.bed"])
        bce.workflow(args)
        result_bed_file_path = "test/test_data/test_bed_file_out.bed"
        result_bed_file = open(result_bed_file_path, 'r')
        interval_list = map(strip, result_bed_file.readlines())

        expected_interval_list = ["chr1:1550-2550",
                                  "chr9:1505-2505",
                                  "chr2:2000-3000"]
        self.assertItemsEqual(expected_interval_list, interval_list)
        os.remove(result_bed_file_path)

    def test_500_length(self):
        interval_param_list = bce.interval_center_extender(self.test_data, 500)

        expected_interval_list = [["chr1", 1800, 2300],
                                  ["chr9", 1755, 2255],
                                  ["chr2", 2250, 2750]]
        self.assertItemsEqual(expected_interval_list, interval_param_list)

    def test_none_length(self):
        interval_param_list = bce.interval_center_extender(self.test_data, None)

        expected_interval_list = [["chr1", 1550, 2550],
                                  ["chr9", 1505, 2505],
                                  ["chr2", 2000, 3000]]
        self.assertItemsEqual(expected_interval_list, interval_param_list)

    def tearDown(self):
        self.test_data.close()
        super(Test, self).tearDown()


if __name__ == "__main__":
    unittest.main()