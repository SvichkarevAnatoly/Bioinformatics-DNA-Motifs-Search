from Bio import motifs
import os
import unittest
import cStringIO
import sys
import MOODS

import lib
import suite
import utils.pattern_matching.pattern_matching as pm


TEST_FASTA_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/fasta.fa")
TEST_PWM_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/pwms_transfac.dat")
TEST_OUT_FILENAME = os.path.join(os.path.dirname(__file__), "test_data/output.dat")


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(Test, cls).setUpClass()
        cls.parser = pm.create_parser()
        cls.seq_search_result = cls.create_seq_search_result()
        cls.pwm_matrix_ctcf = [
            [ 65, 161,  41, 277],  # 1
            [113,  82, 257,  92],  # 2
            [175,  22, 269,  78],  # 3
            [ 32, 481,  14,  17],  # 4
            [  0, 544,   0,   0],  # 5
            [437,   3,  39,  65],  # 6
            [ 17, 304, 216,   7],  # 7
            [ 62, 278,  22, 182],  # 8
            [520,   0,  15,   9],  # 9
            [  0,   0, 544,   0],  # 10
            [220,   3, 318,   3],  # 11
            [ 33,   6, 300, 205],  # 12
            [  5,   0, 536,   3],  # 13
            [ 42,   2, 464,  36],  # 14
            [ 58, 441,   1,  44],  # 15
            [230,   4, 298,  12],  # 16
            [ 47, 298, 175,  24],  # 17
            [ 72, 205,  41, 226],  # 18
            [248,  98, 168,  30]   # 19
        ]

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
        suite.silent_remove(TEST_OUT_FILENAME)
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
            " CAC"
            " 10;17;19;28;35;43;53;61;67"
            " AA",
            "[seq2]"
            " 0;1;6;7;15;19;20;28;31;32;36;44;45;46"
            " TAC"
            " 0;1;2;9;13;34;42"
            " CA",
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
        expected_subseq = "CCC"
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
        expected_subseq = "TTTAAACCCGGGTACG"
        self.assertEqual(len(expected_subseq), len(actual_subseq))
        self.assertEqual(expected_subseq, actual_subseq)

    def test_reverse_complement_excel_best_match_seq(self):
        # reverse complement for GGGGGG
        pwm_str = suite.generate_simple_pwm_str("motif1", "CCCCCC")
        args = suite.create_args("TTGGGGGGACTTGA", pwm_str, True, True, 1.0)

        result = pm.process(args)
        pm.save(result, args)

        actual = suite.read_output_file(args.output)
        expected = "[seq] 2(-) CCCCCC\n"
        self.assertEqual(expected, actual)

    def test_direct_best_match_seq(self):
        pwm_str = suite.generate_simple_pwm_str("motif1", "ACGTAAA")
        args = suite.create_args("ACGTAAA", pwm_str)

        result = pm.process(args)
        pm.save(result, args)

        actual_file_contents = suite.read_output_file(args.output)
        expected_contents = "[seq] 0 ACGTAAA\n"
        self.assertEqual(expected_contents, actual_file_contents)

    def test_revers_complement_best_match_error_score(self):
        sequence = "GAGCGCCACCTGGTGGAGA"
        motif_name = "CTCF"
        pwm_str = suite.generate_pwm_str(motif_name, self.pwm_matrix_ctcf)
        args = suite.create_args(sequence, pwm_str, reverse_complement=True)

        result = pm.process(args)
        pm.save(result, args)

        expected_best_sequence = "CTCGCGGTGGACCACCTCT"
        self.assertEqual(expected_best_sequence, lib.complement(sequence))

        actual_file_contents = suite.read_output_file(args.output)
        expected_contents = "[seq] 0(-) " + expected_best_sequence + "\n"
        self.assertEqual(expected_contents, actual_file_contents)

        pwm = suite.create_pwm(pwm_str)
        matrices = lib.create_matrices_from_pwms(pwm, [motif_name])
        matrix = matrices[0]

        max_score = MOODS.max_score(matrix)
        expected_max_score = 7040
        self.assertEqual(expected_max_score, max_score)
        self.assertEqual(4928, 0.7 * max_score)

        moods_score = result[0].tf_dict[motif_name][0][1]
        self.assertEqual(6429, moods_score)
        self.assertGreaterEqual(moods_score, 0.7 * max_score)

        score = suite.get_score(expected_best_sequence, matrix)
        self.assertEqual(2801, score)
        self.assertGreaterEqual(score, 0.7 * max_score)

    def test_forward_best_match_seq_score_threshold(self):
        sequence = "CTCGCGGTGGACCACCTCT"
        motif_name = "ctcf"
        pwm_str = suite.generate_pwm_str(motif_name, self.pwm_matrix_ctcf)
        args = suite.create_args(sequence, pwm_str)

        result = pm.process(args)
        pm.save(result, args)

        actual_file_contents = suite.read_output_file(args.output)
        expected_best_sequence = "CTCGCGGTGGACCACCTCT"
        expected_contents = "[seq]\n"
        self.assertEqual(expected_contents, actual_file_contents)

        pwm = suite.create_pwm(pwm_str)
        matrices = lib.create_matrices_from_pwms(pwm, [motif_name.upper()])
        matrix = matrices[0]

        max_score = MOODS.max_score(matrix)
        expected_max_score = 7040
        self.assertEqual(expected_max_score, max_score)
        self.assertEqual(4928, 0.7 * max_score)

        score = suite.get_score(expected_best_sequence, matrix)

        self.assertEqual(2801, score)
        self.assertTrue(score <= 0.7 * max_score)

    def assertEqualsMatches(self, args, tf_name, expected_matches):
        result = pm.process(args)[0]

        matches = result.tf_dict[tf_name]
        self.assertItemsEqual(expected_matches, matches)

    def test_direct_MOODS_scores(self):
        pwm_matrix = [
            [100, 200, 300, 400],  # 1
            [ 10,  20,  30,  40],  # 2
            [  1,   2,   3,   4],  # 3
        ]
        tf_name = "CTCF"
        pwm_str = suite.generate_pwm_str(tf_name, pwm_matrix)

        sequence_score_dict = {
            "AAA": [(0, 111)],
            "AAC": [(0, 112)],
            "TGG": [(0, 433)],
        }

        for seq, matches in sequence_score_dict.iteritems():
            args = suite.create_args(seq, pwm_str, threshold=0.0)
            self.assertEqualsMatches(args, tf_name, matches)

    def test_reverse_complement_MOODS_scores(self):
        pwm_matrix = [
            [100, 200, 300, 400],  # 1
            [ 10,  20,  30,  40],  # 2
            [  1,   2,   3,   4],  # 3
        ]
        tf_name = "CTCF"
        pwm_str = suite.generate_pwm_str(tf_name, pwm_matrix)

        sequence_score_dict = {
            "AAA": [(0, 111), (-3, 444)],
        }

        for seq, matches in sequence_score_dict.iteritems():
            args = suite.create_args(seq, pwm_str, reverse_complement=True, threshold=0.0)
            self.assertEqualsMatches(args, tf_name, matches)


if __name__ == "__main__":
    unittest.main()