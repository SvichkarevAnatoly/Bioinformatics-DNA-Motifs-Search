import unittest
import suite


class TestLib(unittest.TestCase):
    def test_nucleotide_to_index(self):
        self.assertEqual(0, suite.to_ind('A'))
        self.assertEqual(1, suite.to_ind('C'))
        self.assertEqual(2, suite.to_ind('G'))
        self.assertEqual(3, suite.to_ind('T'))

    def test_generate_simple_pwm_str(self):
        expected = '\n'.join([
            "ID motif1",
            "P0  A C G T",
            "1   9 0 0 0",
            "2   0 9 0 0",
            "3   0 0 9 0",
            "4   0 0 0 9",
            "//"
        ]) + '\n'
        actual = suite.generate_simple_pwm_str("motif1", "ACGT")
        self.assertEqual(expected, actual)

if __name__ == "__main__":
    unittest.main()
