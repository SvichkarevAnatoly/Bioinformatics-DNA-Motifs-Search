import os
import unittest

import src.bed_center_extender as bce


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), 'test_data/test_bed_file.bed')


class Test(unittest.TestCase):
    def setUp(self):
        super(Test, self).setUp()
        self.test_data = open(TEST_DATA_FILENAME, 'r')

    def test_on_bed_file(self):
        args_list = [self.test_data, 500, None]
        bed_interval_param_list = bce.bed_center_extender(args_list)

        expected_interval = ["chr1", 1800, 2300]
        self.assertEqual(expected_interval, bed_interval_param_list[0])

    def tearDown(self):
        self.test_data.close()
        super(Test, self).tearDown()

if __name__ == "__main__":
    unittest.main()