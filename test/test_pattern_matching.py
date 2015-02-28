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

    def test_reversed_flag(self):
        args = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)]
        self.args = self.parser.parse_args(args)
        self.assertFalse(self.args.reversed)

        args.append("--reversed")
        self.args = self.parser.parse_args(args)
        self.assertTrue(self.args.reversed)

    def test_run_on_test_data(self):
        self.args = self.parser.parse_args([str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)])
        result = pm.process(self.args)

        self.assertEqual(2, len(result))
        seq1_result = result[0]
        self.assertEqual("seq1", seq1_result[0])
        seq1_tf1_result = seq1_result[1][0]
        self.assertEqual("motif1", seq1_tf1_result[0])
        expected_pos = (0, 6.0)
        seq1_tf1_positions = seq1_tf1_result[1]
        self.assertEqual(expected_pos, seq1_tf1_positions[0])

    def tearDown(self):
        super(Test, self).tearDown()

    @classmethod
    def tearDownClass(cls):
        super(Test, cls).tearDownClass()


if __name__ == "__main__":
    unittest.main()