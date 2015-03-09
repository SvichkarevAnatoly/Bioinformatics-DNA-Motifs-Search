import os
import unittest
import errno
import cStringIO

import utils.pattern_matching as pm


TEST_FASTA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/fasta.fa")
TEST_PWM_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/pwms_transfac.dat")
TEST_OUT_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/output.dat")


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

    def test_reverse_complement_flag(self):
        args = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)]
        self.args = self.parser.parse_args(args)
        self.assertFalse(self.args.reverse_complement)

        args.append("--reverse-complement")
        self.args = self.parser.parse_args(args)
        self.assertTrue(self.args.reverse_complement)

        tempfile = cStringIO.StringIO()
        self.args.output = tempfile

        result = pm.process(self.args)
        pm.save(result, self.args)

        tempfile.seek(0)
        actual_file_contents = tempfile.read()
        expected_contents = "\n".join([
            ">seq1",
            "motif1 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65;"
            "-66;-63;-57;-56;-49;-48;-45;-39;-38;-35;-34;-33;-27;-20;-14;-3",
            "motif2 10;17;19;28;35;43;53;61;67;-69;-66;-57;-50",
            ">seq2",
            "motif1 0;1;6;7;15;19;20;28;31;32;36;44;45;46;"
            "-46;-39;-38;-32;-31;-27;-25;-14;-10;-9;-6",
            "motif2 0;1;2;9;13;34;42;-44;-32;-31;-19;-14;-11",
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

    def test_backward_flag(self):
        args = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)]
        self.args = self.parser.parse_args(args)
        self.assertFalse(self.args.backward)

        args.append("--backward")
        self.args = self.parser.parse_args(args)
        self.assertTrue(self.args.backward)

        tempfile = cStringIO.StringIO()
        self.args.output = tempfile

        result = pm.process(self.args)
        pm.save(result, self.args)

        tempfile.seek(0)
        actual_file_contents = tempfile.read()
        expected_contents = "\n".join([
            ">seq1",
            "motif1 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65",
            "backward 2;3;8;9;12;16;19;20;26;27;34;41;48;49;50;52;55;56;59;64;67;68",
            "motif2 10;17;19;28;35;43;53;61;67",
            "backward 16;35;41;50;52;53",
            ">seq2",
            "motif1 0;1;6;7;15;19;20;28;31;32;36;44;45;46",
            "backward 0;1;5;9;12;13;14;17;18;26;30;31;34;38;39;45;46;47",
            "motif2 0;1;2;9;13;34;42",
            "backward 32;46;47"
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

    def test_backward_excel_flag(self):
        args_str = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME),
                    "--backward", "--excel", "-tf", "motif1"]
        args = self.parser.parse_args(args_str)
        self.assertTrue(args.backward)
        self.assertTrue(args.backward)

        tempfile = cStringIO.StringIO()
        args.output = tempfile

        result = pm.process(args)
        pm.save(result, args)

        tempfile.seek(0)
        actual_file_contents = tempfile.read()
        expected_contents = "\n".join([
            "[seq1] "
            "0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65|"
            "2;3;8;9;12;16;19;20;26;27;34;41;48;49;50;52;55;56;59;64;67;68",
            "[seq2] "
            "0;1;6;7;15;19;20;28;31;32;36;44;45;46|"
            "0;1;5;9;12;13;14;17;18;26;30;31;34;38;39;45;46;47"
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

    def test_run_on_test_data(self):
        self.args = self.parser.parse_args([str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)])
        result = pm.process(self.args)

        expected_seq_number = 2
        self.assertEqual(expected_seq_number, len(result))
        seq_result = result[0]
        self.assertEqual("seq1", seq_result.seq_name)
        expected_tf = "motif1"
        self.assertEqual(expected_tf, seq_result.tfs[0])
        self.assertIsNotNone(seq_result.tf_dict[expected_tf])

        expected_first_match = (0, 6.0)
        actual_matching_list = seq_result.tf_dict[expected_tf].directed
        self.assertEqual(expected_first_match, actual_matching_list[0])

    def test_several_tfs(self):
        args = self.parser.parse_args([
            str(TEST_FASTA_FILENAME),
            str(TEST_PWM_FILENAME),
            "-tf",
            "motif1"
        ])
        self.assertEqual(1, len(args.tf))
        self.assertEqual(["motif1"], args.tf)

        result = pm.process(args)

        expected_seq_number = 2
        self.assertEqual(expected_seq_number, len(result))
        seq_result = result[0]
        self.assertEqual("seq1", seq_result.seq_name)
        expected_tf = "motif1"
        self.assertEqual(expected_tf, seq_result.tfs[0])
        self.assertIsNotNone(seq_result.tf_dict[expected_tf])

        expected_first_match = (0, 6.0)
        actual_matching_list = seq_result.tf_dict[expected_tf].directed
        self.assertEqual(expected_first_match, actual_matching_list[0])

    def test_output_file_result(self):
        self.args = self.parser.parse_args([str(TEST_FASTA_FILENAME),
                                            str(TEST_PWM_FILENAME),
                                            "-o",
                                            str(TEST_OUT_FILENAME)])
        self.assertEquals(TEST_OUT_FILENAME, self.args.output.name)
        silent_remove(TEST_OUT_FILENAME)
        tempfile = cStringIO.StringIO()
        self.args.output = tempfile

        result = pm.process(self.args)
        pm.save(result, self.args)

        tempfile.seek(0)
        actual_file_contents = tempfile.read()

        expected_contents = '\n'.join([
            ">seq1",
            "motif1 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65",
            "motif2 10;17;19;28;35;43;53;61;67",
            ">seq2",
            "motif1 0;1;6;7;15;19;20;28;31;32;36;44;45;46",
            "motif2 0;1;2;9;13;34;42",
        ]) + '\n'

        self.assertEqual(expected_contents, actual_file_contents)

    def test_save_excel_format(self):
        self.args = self.parser.parse_args([str(TEST_FASTA_FILENAME),
                                            str(TEST_PWM_FILENAME),
                                            "-e"])
        tempfile = cStringIO.StringIO()
        self.args.output = tempfile

        result = pm.process(self.args)
        pm.save(result, self.args)

        tempfile.seek(0)
        actual_file_contents = tempfile.read()

        expected_contents = '\n'.join([
            "[seq1] 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65 10;17;19;28;35;43;53;61;67",
            "[seq2] 0;1;6;7;15;19;20;28;31;32;36;44;45;46 0;1;2;9;13;34;42",
        ]) + '\n'

        self.assertEqual(expected_contents, actual_file_contents)

    def tearDown(self):
        super(Test, self).tearDown()

    @classmethod
    def tearDownClass(cls):
        super(Test, cls).tearDownClass()


if __name__ == "__main__":
    unittest.main()