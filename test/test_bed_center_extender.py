import os
from strop import strip
import unittest
import errno

import src.bed_center_extender as bce


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/bed_file.bed")
TEST_DATA_OUT_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/bed_file_out.bed")


def silent_remove(file_name):
    try:
        os.remove(file_name)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


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
        args = self.parser.parse_args([str(TEST_DATA_FILENAME)])
        bce.workflow(args)
        result_bed_file = open(TEST_DATA_OUT_FILENAME, 'r')
        interval_list = map(strip, result_bed_file.readlines())

        expected_interval_list = ["chr1:1550-2550",
                                  "chr9:1505-2505",
                                  "chr2:2000-3000"]
        self.assertItemsEqual(expected_interval_list, interval_list)

    def test_only_bedfile_and_500_length_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME), "-l500"])
        bce.workflow(args)
        result_bed_file = open(TEST_DATA_OUT_FILENAME, 'r')
        interval_list = map(strip, result_bed_file.readlines())

        expected_interval_list = ["chr1:1800-2300",
                                  "chr9:1755-2255",
                                  "chr2:2250-2750"]
        self.assertItemsEqual(expected_interval_list, interval_list)

    def test_only_bedfile_and_outfile_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME), "-o", str(TEST_DATA_OUT_FILENAME)])
        bce.workflow(args)
        result_bed_file = open(TEST_DATA_OUT_FILENAME, 'r')
        interval_list = map(strip, result_bed_file.readlines())

        expected_interval_list = ["chr1:1550-2550",
                                  "chr9:1505-2505",
                                  "chr2:2000-3000"]
        self.assertItemsEqual(expected_interval_list, interval_list)

    def test_full_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME), "-l500", "-o", str(TEST_DATA_OUT_FILENAME)])
        bce.workflow(args)
        result_bed_file = open(TEST_DATA_OUT_FILENAME, 'r')
        interval_list = map(strip, result_bed_file.readlines())

        expected_interval_list = ["chr1:1800-2300",
                                  "chr9:1755-2255",
                                  "chr2:2250-2750"]
        self.assertItemsEqual(expected_interval_list, interval_list)

    def test_wrong_cl_args_and_not_creating_file(self):
        silent_remove(TEST_DATA_OUT_FILENAME)

        args_list = [str(TEST_DATA_FILENAME), "-l", "not_number", "-o", str(TEST_DATA_OUT_FILENAME)]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args_list)

        self.assertTrue(not os.path.isfile(TEST_DATA_OUT_FILENAME))

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

    @classmethod
    def tearDownClass(cls):
        silent_remove(TEST_DATA_OUT_FILENAME)
        super(Test, cls).tearDownClass()


if __name__ == "__main__":
    unittest.main()