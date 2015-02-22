import os
import unittest

import src.bed_center_extender as bce


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), 'test_data/test_bed_file.bed')


class Test(unittest.TestCase):
    def setUp(self):
        self.test_data = open(TEST_DATA_FILENAME, 'r')

    def test_500_length(self):
        args_list = [self.test_data, 500, None]
        bed_interval_param_list = bce.bed_center_extender(args_list)

        expected_interval_list = [["chr1", 1800, 2300],
                                  ['chr9', 1755, 2255],
                                  ['chr2', 2250, 2750]]
        self.assertItemsEqual(expected_interval_list, bed_interval_param_list)

    def test_none_length(self):
        args_list = [self.test_data, None, None]
        bed_interval_param_list = bce.bed_center_extender(args_list)

        expected_interval_list = [["chr1", 1550, 2550],
                                  ['chr9', 1505, 2505],
                                  ['chr2', 2000, 3000]]
        self.assertItemsEqual(expected_interval_list, bed_interval_param_list)

    def tearDown(self):
        self.test_data.close()
        super(Test, self).tearDown()

if __name__ == "__main__":
    unittest.main()