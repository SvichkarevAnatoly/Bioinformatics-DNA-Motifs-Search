import unittest

import analysis.excel_ucsc_ids_to_bed_order as eu


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = eu.create_parser()

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.args = self.parser.parse_args([])


if __name__ == "__main__":
    unittest.main()