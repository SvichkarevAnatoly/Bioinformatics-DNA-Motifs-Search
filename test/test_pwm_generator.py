from argparse import Namespace
import unittest
import cStringIO

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import utils.pwm_generator as pg


class Test(unittest.TestCase):
    def create_seqs(self, sequences):
        return [Seq(seq, IUPAC.unambiguous_dna) for seq in sequences]

    def create_args(self, sequences):
        args = Namespace()
        args.seqs = self.create_seqs(sequences)
        args.output = cStringIO.StringIO()
        return args

    def test_creating_motif(self):
        seqs = [
            "ACGT",
            "AAAA"
        ]
        args = self.create_args(seqs)
        motif = pg.process(args)

        motif_seqs = map(str, motif.instances)
        self.assertItemsEqual(seqs, motif_seqs)

    def test_check_raise_value_error_diff_len(self):
        seqs = [
            "ACGT",
            "A"
        ]
        args = self.create_args(seqs)

        with self.assertRaisesRegexp(ValueError, "Seqs length differ"):
            pg.check(args.seqs)

    def test_check_raise_value_error_alphabet(self):
        seqs = [
            "ACGX",
        ]
        args = self.create_args(seqs)

        with self.assertRaisesRegexp(ValueError, "Seq contains letter not from Alphabet"):
            pg.check(args.seqs)

if __name__ == "__main__":
    unittest.main()
