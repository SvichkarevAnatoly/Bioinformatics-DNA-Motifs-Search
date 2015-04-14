import os
import unittest

from bbcflib.track import track


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), "../data_samples/analysis/csv_like.bed")


class Test(unittest.TestCase):
    def test_track_fields(self):
        t = track(TEST_DATA_FILENAME, format='txt', separator='\t',
                  fields=['chrom', 'chromStart', 'chromEnd',
                          'c1', 'c2', 'c3', 'c4', 'c5', 'c6',
                          'peakName'])

        s = t.read(fields=['chrom', 'chromStart', 'chromEnd', 'peakName'])

        peak1 = s.next()
        expected_peak1 = ("chr1", "4736010", "4736158", "Z4_Sox2_peak_1")
        self.assertSequenceEqual(expected_peak1, peak1)

        peak2 = s.next()
        expected_peak2 = ("chr1", "5223047", "5223196", "Z4_Sox2_peak_2")
        self.assertSequenceEqual(expected_peak2, peak2)

if __name__ == "__main__":
    unittest.main()