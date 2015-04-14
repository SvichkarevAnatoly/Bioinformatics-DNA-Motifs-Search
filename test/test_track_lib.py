import os
import unittest

from bbcflib.track import track


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bed_file_path = os.path.join(os.path.dirname(__file__),
                                         "../data_samples/analysis/csv_like.bed")
        cls.bed_fields = ['chrom', 'chromStart', 'chromEnd',
                          'c1', 'c2', 'c3', 'c4', 'c5', 'c6',
                          'peakName']

    def test_track_fields(self):
        t = track(self.bed_file_path, format='txt',
                  separator='\t', fields=self.bed_fields)

        s = t.read(fields=['chrom', 'chromStart', 'chromEnd', 'peakName'])

        peak1 = s.next()
        expected_peak1 = ("chr1", "4736010", "4736158", "Z4_Sox2_peak_1")
        self.assertSequenceEqual(expected_peak1, peak1)

        peak2 = s.next()
        expected_peak2 = ("chr1", "5223047", "5223196", "Z4_Sox2_peak_2")
        self.assertSequenceEqual(expected_peak2, peak2)


if __name__ == "__main__":
    unittest.main()