import os
import unittest

from bbcflib.track import track
from bbcflib.track import FeatureStream
import suite


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bed_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/csv_like.bed")
        cls.bed_out_file_path = os.path.join(
            os.path.dirname(__file__),
            "../data_samples/analysis/csv_like_out.bed")

        cls.bed_fields = ['chrom', 'chromStart', 'chromEnd',
                          'c1', 'c2', 'c3', 'c4', 'c5', 'c6',
                          'peakName']

    def tearDown(self):
        suite.silent_remove(self.bed_out_file_path)

    def test_read_csv_fields(self):
        t = track(self.bed_file_path, format='txt',
                  separator='\t', fields=self.bed_fields)

        s = t.read(fields=['chrom', 'chromStart', 'chromEnd', 'peakName'])

        expected_peak1 = ("chr1", "4736010", "4736158", "Z4_Sox2_peak_1")
        self.assertSequenceEqual(expected_peak1, s.next())

        expected_peak2 = ("chr1", "5223047", "5223196", "Z4_Sox2_peak_2")
        self.assertSequenceEqual(expected_peak2, s.next())

    def test_simple_write_bed_fields(self):
        peaks = [
            ('chr1', 12, 13),
            ('chr1', 24, 28)]

        s_write = FeatureStream(peaks, fields=['chr', 'start', 'end'])
        t_write = track(self.bed_out_file_path, fields=s_write.fields)
        t_write.write(s_write)

        t_read = track(self.bed_out_file_path)
        s_read = t_read.read()

        self.assertSequenceEqual(peaks[0], s_read.next())
        self.assertSequenceEqual(peaks[1], s_read.next())


if __name__ == "__main__":
    unittest.main()