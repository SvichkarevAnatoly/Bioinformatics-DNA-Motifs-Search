from argparse import Namespace
from Bio import motifs
from Bio import SeqIO
import os
import unittest
import errno
import cStringIO
import sys

import utils.pattern_matching.pattern_matching as pm


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
        cls.seq_search_result = cls.create_seq_search_result()

    @classmethod
    def create_seq_search_result(cls):
        tf_name = "NANOG"
        sequence = "AAATTTGGGCCCATGC"
        # complem  "TTTAAACCCGGGTACG"
        # rev_com  "GCATGGGCCCAAATTT"
        seq_search_result = pm.SeqSearchResults("", sequence, [tf_name], [3])
        return seq_search_result

    def test_with_empty_args(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_default_args(self):
        args = self.parser.parse_args([str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)])
        self.assertEqual(sys.stdout, args.output)
        self.assertEqual(None, args.tf)
        self.assertEqual(0.7, args.threshold)
        self.assertFalse(args.reverse_complement)
        self.assertFalse(args.excel)

    def test_reverse_complement_flag(self):
        args_str = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)]
        args = self.parser.parse_args(args_str)
        self.assertFalse(args.reverse_complement)

        args_str.append("--reverse-complement")
        args = self.parser.parse_args(args_str)
        self.assertTrue(args.reverse_complement)

        args.output = cStringIO.StringIO()

        result = pm.process(args)
        pm.save(result, args)

        args.output.seek(0)
        actual_file_contents = args.output.read()
        expected_contents = "\n".join([
            ">seq1",
            "MOTIF1 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;"
            "62;65;6(-);9(-);15(-);16(-);23(-);24(-);27(-);33(-);"
            "34(-);37(-);38(-);39(-);45(-);52(-);58(-);69(-)",
            "MOTIF2 10;17;19;28;35;43;53;61;67;3(-);6(-);15(-);22(-)",
            ">seq2",
            "MOTIF1 0;1;6;7;15;19;20;28;31;32;36;44;45;46;4(-);11(-);"
            "12(-);18(-);19(-);23(-);25(-);36(-);40(-);41(-);44(-)",
            "MOTIF2 0;1;2;9;13;34;42;6(-);18(-);19(-);31(-);36(-);39(-)",
        ]) + '\n'
        self.assertEqual(expected_contents, actual_file_contents)

    def test_run_on_test_data(self):
        args = self.parser.parse_args([str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME)])
        result = pm.process(args)

        expected_seq_number = 2
        self.assertEqual(expected_seq_number, len(result))
        seq_result = result[0]
        self.assertEqual("seq1", seq_result.seq_name)
        expected_tf = "MOTIF1"
        self.assertEqual(expected_tf, seq_result.tfs[0])
        self.assertIsNotNone(seq_result.tf_dict[expected_tf])

        expected_first_match = (0, 6.0)
        actual_matches_list = seq_result.tf_dict[expected_tf]
        self.assertEqual(expected_first_match, actual_matches_list[0])

    def test_upper_case_tfs(self):
        args_str1 = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME), "-tf", "motif1"]
        args1 = self.parser.parse_args(args_str1)
        self.assertEqual(["MOTIF1"], args1.tf)

        args_str2 = [str(TEST_FASTA_FILENAME), str(TEST_PWM_FILENAME), "-tf", "MOTIF1"]
        args2 = self.parser.parse_args(args_str2)
        self.assertEqual(["MOTIF1"], args2.tf)

    def test_several_tfs(self):
        args = self.parser.parse_args([
            str(TEST_FASTA_FILENAME),
            str(TEST_PWM_FILENAME),
            "-tf",
            "motif1"
        ])
        self.assertEqual(1, len(args.tf))
        self.assertEqual(["MOTIF1"], args.tf)

        result = pm.process(args)

        expected_seq_number = 2
        self.assertEqual(expected_seq_number, len(result))
        seq_result = result[0]
        self.assertEqual("seq1", seq_result.seq_name)
        expected_tf = "MOTIF1"
        self.assertEqual(expected_tf, seq_result.tfs[0])
        self.assertIsNotNone(seq_result.tf_dict[expected_tf])

        expected_first_match = (0, 6.0)
        actual_matches_list = seq_result.tf_dict[expected_tf]
        self.assertEqual(expected_first_match, actual_matches_list[0])

    def test_output_file_result(self):
        args = self.parser.parse_args([str(TEST_FASTA_FILENAME),
                                       str(TEST_PWM_FILENAME),
                                       "-o",
                                       str(TEST_OUT_FILENAME)])
        self.assertEquals(TEST_OUT_FILENAME, args.output.name)
        silent_remove(TEST_OUT_FILENAME)
        args.output = cStringIO.StringIO()

        result = pm.process(args)
        pm.save(result, args)

        args.output.seek(0)
        actual_file_contents = args.output.read()

        expected_contents = '\n'.join([
            ">seq1",
            "MOTIF1 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65",
            "MOTIF2 10;17;19;28;35;43;53;61;67",
            ">seq2",
            "MOTIF1 0;1;6;7;15;19;20;28;31;32;36;44;45;46",
            "MOTIF2 0;1;2;9;13;34;42",
        ]) + '\n'

        self.assertEqual(expected_contents, actual_file_contents)

    def test_save_excel_format(self):
        args = self.parser.parse_args([str(TEST_FASTA_FILENAME),
                                       str(TEST_PWM_FILENAME),
                                       "-e"])
        args.output = cStringIO.StringIO()

        result = pm.process(args)
        pm.save(result, args)

        args.output.seek(0)
        actual_file_contents = args.output.read()

        expected_contents = '\n'.join([
            "[seq1]"
            " 0;3;6;11;12;15;16;17;19;28;34;41;48;49;53;58;59;62;65"
            " GGGCGCACGCGGA"
            " 10;17;19;28;35;43;53;61;67"
            " CGCGGAAGGGAG",
            "[seq2]"
            " 0;1;6;7;15;19;20;28;31;32;36;44;45;46"
            " AGGCATACTTTCC"
            " 0;1;2;9;13;34;42"
            " CCTGCCATTCTG",
        ]) + '\n'

        self.assertEqual(expected_contents, actual_file_contents)

    def test_seq_search_results_best_match(self):
        tf_name = "NANOG"
        seq_search_result = pm.SeqSearchResults("", "", [tf_name], [3])
        matches = [(3, 4.0), (1, 6.0), (-6, 3.0)]
        seq_search_result.fill_matches([matches])

        actual_best_match = seq_search_result.best_match(tf_name)
        expected_best_match = matches[1]
        self.assertEqual(expected_best_match, actual_best_match)

    def test_seq_search_results_nearest_to_center_best_match(self):
        tf_name = "NANOG"
        sequence = "AAAAAAAAAA"
        seq_search_result = pm.SeqSearchResults("", sequence, [tf_name], [3])
        matches = [(3, 4.0), (1, 6.0), (-6, 3.0), (4, 6.0)]
        seq_search_result.fill_matches([matches])

        actual_best_match = seq_search_result.best_match(tf_name)
        expected_best_match = matches[3]
        self.assertEqual(expected_best_match, actual_best_match)

    def test_seq_search_results_match_subseq_positive_zero_delta(self):
        tf_name = "NANOG"
        match = (6, 2.0)
        delta = 0
        actual_subseq = self.seq_search_result.match_subseq(match[0], tf_name, delta)
        expected_subseq = "GGG"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def test_seq_search_results_match_subseq_negative_zero_delta(self):
        tf_name = "NANOG"
        match = (-10, 2.0)  # -10 + 16 = 6
        delta = 0
        actual_subseq = self.seq_search_result.match_subseq(match[0], tf_name, delta)
        expected_subseq = "GCC"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def test_seq_search_results_match_subseq_positive_one_delta(self):
        tf_name = "NANOG"
        match = (6, 2.0)
        delta = 1
        actual_subseq = self.seq_search_result.match_subseq(match[0], tf_name, delta)
        expected_subseq = "TGGGC"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def test_seq_search_results_match_subseq_positive_end_range(self):
        tf_name = "NANOG"
        match = (6, 2.0)
        delta = 10
        actual_subseq = self.seq_search_result.match_subseq(match[0], tf_name, delta)
        expected_subseq = "AAATTTGGGCCCATGC"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def test_seq_search_results_match_subseq_negative_end_range(self):
        tf_name = "NANOG"
        match = (-10, 2.0)
        delta = 10
        actual_subseq = self.seq_search_result.match_subseq(match[0], tf_name, delta)
        expected_subseq = "GCATGGGCCCAAATTT"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def fill_args(self):
        args = Namespace()
        args.pwm = self.create_pwm()
        args.fasta = self.create_fasta()
        args.output = cStringIO.StringIO()
        args.tf = ["MOTIF1"]
        args.reverse_complement = True
        args.excel = True
        args.threshold = 0.7

        return args

    def read_output_file(self, output):
        output.seek(0)
        actual_file_contents = output.read()
        return actual_file_contents

    def test_reverse_complement_excel_best_match_seq(self):
        args = self.fill_args()

        result = pm.process(args)
        pm.save(result, args)

        actual_file_contents = self.read_output_file(args.output)

        expected_contents = "[seq2]" \
                            " 0;1;6;7;15;19;20;28;31;32;36;44;45;46" \
                            ";4(-);11(-);12(-);18(-);19(-);23(-);25(-);36(-);40(-);41(-);44(-)" \
                            " AGGCATACTTTCC\n"

        self.assertEqual(expected_contents, actual_file_contents)

    def create_pwm(self):
        pwm_str = '\n'.join([
            "ID  motif1",
            "P0      A      C      G      T",
            "01      1      1      1      2      N",
            "02      2      1      0      1      A",
            "03      3      4      0      0      M",
            "//"
        ]) + '\n'
        pwm_handler = cStringIO.StringIO()
        pwm_handler.write(pwm_str)
        pwm_handler.seek(0)
        pwm_records = motifs.parse(pwm_handler, "TRANSFAC")
        pwm_handler.close()
        return pwm_records

    def create_fasta(self):
        fasta_str = '\n'.join([
            ">seq2",
            "CAAAGATGCCAGGCATACTTTCCGAGTAGCCTGCCATTCTGGCAGTCCCG"
        ]) + '\n'
        fasta_handler = cStringIO.StringIO()
        fasta_handler.write(fasta_str)
        fasta_handler.seek(0)
        fasta = list(SeqIO.parse(fasta_handler, "fasta"))
        fasta_handler.close()
        return fasta


if __name__ == "__main__":
    unittest.main()