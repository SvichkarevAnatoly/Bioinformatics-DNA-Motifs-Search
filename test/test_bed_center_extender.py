import os
import unittest
import cStringIO

import utils.bed_center_extender as bce


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/bed_file.bed")


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = bce.create_parser()
        cls.input_intervals = [
            "chr1:2000-2100",
            "chr9:2000-2010",
            "chr2:2000-3000"
        ]

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_only_bedfile_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME)])

        result = bce.process(args)
        expected_interval_list = [
            ["chr1", 1550, 2550],
            ["chr9", 1505, 2505],
            ["chr2", 2000, 3000]
        ]
        self.assertItemsEqual(expected_interval_list, result)

    def test_only_bedfile_and_500_length_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME), "-l500"])
        self.assertEquals(500, args.length)

        result = bce.process(args)
        expected_interval_list = [
            ["chr1", 1800, 2300],
            ["chr9", 1755, 2255],
            ["chr2", 2250, 2750]
        ]
        self.assertItemsEqual(expected_interval_list, result)

    def test_only_bedfile_and_outfile_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME)])

        args.output = cStringIO.StringIO()

        result = bce.process(args)
        bce.save(result, args)

        args.output.seek(0)
        actual_contents = args.output.read()
        expected_contents = "\n".join([
            "chr1:1550-2550",
            "chr9:1505-2505",
            "chr2:2000-3000"
        ]) + '\n'
        self.assertItemsEqual(expected_contents, actual_contents)

    def test_full_cl_args(self):
        args = self.parser.parse_args([str(TEST_DATA_FILENAME), "-l500"])
        self.assertEquals(500, args.length)

        args.output = cStringIO.StringIO()

        result = bce.process(args)
        bce.save(result, args)

        args.output.seek(0)
        actual_contents = args.output.read()
        expected_contents = "\n".join([
            "chr1:1800-2300",
            "chr9:1755-2255",
            "chr2:2250-2750"
        ]) + '\n'
        self.assertItemsEqual(expected_contents, actual_contents)

    def test_wrong_cl_args_and_not_creating_file(self):
        args_list = [str(TEST_DATA_FILENAME), "-l", "not_number"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args_list)

    def test_500_length(self):
        interval_param_list = bce.interval_center_extender(self.input_intervals, 500)
        expected_interval_list = [
            ["chr1", 1800, 2300],
            ["chr9", 1755, 2255],
            ["chr2", 2250, 2750]
        ]
        self.assertItemsEqual(expected_interval_list, interval_param_list)

    def test_none_length(self):
        interval_param_list = bce.interval_center_extender(self.input_intervals, None)
        expected_interval_list = [
            ["chr1", 1550, 2550],
            ["chr9", 1505, 2505],
            ["chr2", 2000, 3000]
        ]
        self.assertItemsEqual(expected_interval_list, interval_param_list)

if __name__ == "__main__":
    unittest.main()