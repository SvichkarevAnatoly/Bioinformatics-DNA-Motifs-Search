import sys

if len(sys.argv) == 4:
    input_bed_file_name = sys.argv[1]
    extended_length = sys.argv[2]
    output_bed_file_name = sys.argv[3]
else:
    print "usage: fasta_motif_center.py input_file [length] [output_file]\n" \
          "input_file  : file with bed format intervals\n" \
          "length      : common extended length\n" \
          "output_file : output file with extended bed format intervals"

# TODO: read bed file with intervals different length


# TODO: compute central range 500 length


# TODO: write to output file