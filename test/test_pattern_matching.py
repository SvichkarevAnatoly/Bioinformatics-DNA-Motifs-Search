import os
import unittest
import errno

import src.pattern_matching as pm

TEST_FASTA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/fasta.fa")
TEST_PWM_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/pwms_transfac.dat")


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
            self.args = self.parser.parse_args([])

    def test_positional_args(self):
        self.args = self.parser.parse_args([str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)])
        # TODO: how to check open file
        self.assertEqual(0.7, self.args.threshold)
        # TODO: how to close in tearDown?
        self.args.fasta.close()
        self.args.pwm.close()

    def tearDown(self):
        super(Test, self).tearDown()

    @classmethod
    def tearDownClass(cls):
        super(Test, cls).tearDownClass()


if __name__ == "__main__":
    unittest.main()