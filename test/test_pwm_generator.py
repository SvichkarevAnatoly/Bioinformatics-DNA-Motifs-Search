from argparse import Namespace
import unittest
import cStringIO
import lib

import suite
import utils.pwm_generator as pg


class Test(unittest.TestCase):
    def create_args(self, sequences, motif_name):
        args = Namespace()
        args.seqs = lib.biopython_seqs(sequences)
        args.output = cStringIO.StringIO()
        args.motif = motif_name
        return args

    def test_creating_motif(self):
        seqs = [
            "ACGT",
            "AAAA"
        ]
        args = self.create_args(seqs, None)
        motif = pg.process(args)

        motif_seqs = map(str, motif.instances)
        self.assertItemsEqual(seqs, motif_seqs)

    def test_save(self):
        seqs = [
            "ACGT",
            "AAAA"
        ]
        args = self.create_args(seqs, None)
        motif = pg.process(args)
        pg.save(motif, args)

        expected_file_contents = '\n'.join([
            "P0      A      C      G      T",
            "01      2      0      0      0      A",
            "02      1      1      0      0      M",
            "03      1      0      1      0      R",
            "04      1      0      0      1      W",
            "XX",
            "//"
        ]) + '\n'
        actual_file_contents = suite.read_output_file(args.output)
        self.assertEqual(expected_file_contents, actual_file_contents)

    def test_save_with_name(self):
        seqs = [
            "ACGT",
            "AAAA"
        ]
        args = self.create_args(seqs, "motif1")
        motif = pg.process(args)
        pg.save(motif, args)

        expected_file_contents = '\n'.join([
            "ID  motif1",
            "P0      A      C      G      T",
            "01      2      0      0      0      A",
            "02      1      1      0      0      M",
            "03      1      0      1      0      R",
            "04      1      0      0      1      W",
            "XX",
            "//"
        ]) + '\n'
        actual_file_contents = suite.read_output_file(args.output)
        self.assertEqual(expected_file_contents, actual_file_contents)


if __name__ == "__main__":
    unittest.main()
