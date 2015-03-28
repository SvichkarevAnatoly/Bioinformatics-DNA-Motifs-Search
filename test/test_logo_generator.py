import unittest

import utils.logo_generator as lg


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = lg.create_parser()

    def test_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

if __name__ == "__main__":
    unittest.main()

