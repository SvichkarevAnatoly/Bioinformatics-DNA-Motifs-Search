from argparse import Namespace
import unittest
import cStringIO


class Test(unittest.TestCase):
    def create_seqs_handler(self, sequences):
        seqs_handler = cStringIO.StringIO()
        seqs_handler.write('\n'.join(sequences)+'\n')
        seqs_handler.close()
        return seqs_handler

    def create_args(self, sequences):
        args = Namespace()
        args.seqs = self.create_seqs_handler(sequences)
        args.output = cStringIO.StringIO()
        return args

    def test_creating_pwm(self):
        pass

if __name__ == "__main__":
    unittest.main()
