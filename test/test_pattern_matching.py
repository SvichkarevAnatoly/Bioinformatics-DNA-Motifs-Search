import os
import unittest
import errno

import src.pattern_matching as pm

TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/fasta.fa")
TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/fasta.fa")


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
        cls.parser = pm.create_parser()

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    @classmethod
    def tearDownClass(cls):
        super(Test, cls).tearDownClass()


if __name__ == "__main__":
    unittest.main()