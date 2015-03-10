import argparse
import unittest
import cStringIO

import analysis.excel_ucsc_ids_to_bed_order as eu


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = eu.create_parser()

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def createExcelFile(self, args):
        args.excel = cStringIO.StringIO()
        excel_input = "\n".join([
            "[mm10_ct_UserTrack_3545_0"
            " range=chr1:184025039-184025538"
            " 5'pad=0 3'pad=0 strand=+ repeatMasking=none]"
            " 227;-190;-65|297 #|#",
            "[mm10_ct_UserTrack_3545_1"
            " range=chr1:3062703-3063202"
            " 5'pad=0 3'pad=0 strand=+ repeatMasking=none]"
            " 226;476;-250|-266;-253 34|#",
        ]) + '\n'
        args.excel.write(excel_input)
        args.excel.seek(0)

    def createOutputFile(self, args):
        args.output = cStringIO.StringIO()

    def test_workflow(self):
        args = argparse.Namespace()
        self.createExcelFile(args)
        self.createOutputFile(args)

        result = eu.process(args)
        eu.save(result, args)

        args.output.seek(0)
        expected_file_contents = "\n".join([
            "chr1:184025039-184025538"
            " 227;-190;-65|297 #|#",
            "chr1:3062703-3063202"
            " 226;476;-250|-266;-253 34|#",
        ]) + '\n'
        actual_file_contents = args.output.read()
        self.assertEqual(expected_file_contents, actual_file_contents)


if __name__ == "__main__":
    unittest.main()