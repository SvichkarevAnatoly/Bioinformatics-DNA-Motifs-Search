# Bioinformatics-DNA-Motifs-Search
## Requirements
+ Recommend using [Python 2.7](http://www.python.org).
+ [Biopython](http://biopython.org/)
+ [MOOD](http://www.cs.helsinki.fi/group/pssmfind/)

## Usage
#### bed_center_extender.py [-h] [-l LENGTH] [-o OUTFILE] bedfile

    Resulting intervals in the input file to the same length

    positional arguments:
    bedfile               file with bed format intervals

    optional arguments:
    -h, --help            show this help message and exit
    -l LENGTH, --length LENGTH
                        common extended length. If not specified, is extended
                        to the maximum length of the interval in the input
                        file.
    -o OUTFILE, --output OUTFILE
                        output file with extended bed format intervals

#### pattern_matching.py [-h] [-o [OUTPUT]] [-tf TF [TF ...]] [-th THRESHOLD] [-rc] [-b] [-e] fasta pwm
    
    Matching position weight matrices (PWM) against DNA sequences
    
    positional arguments:
      fasta                 fasta file with DNA sequences
      pwm                   file with position weight matrices (PWM)
    
    optional arguments:
      -h, --help            show this help message and exit
      -o [OUTPUT], --output [OUTPUT]
                            output file with matching results. Default stdout.
      -tf TF [TF ...], --factor TF [TF ...]
                            transcription factor name in pwm file. Default
                            matching with all tf in pwm file.
      -th THRESHOLD, --threshold THRESHOLD
                            The parameter threshold split for better control on
                            what parts of the scoring are used. Default 0.7.
      -rc, --reverse-complement
                            Scans against reverse complement sequence in addition
                            to the input sequence. Hits on reverse complement are
                            reported at position [position - sequence_length],
                            which is always negative. The actual hit site for any
                            hit is always seq[pos, pos + matrix_length]. Default
                            False.
      -b, --backward        For searching in both direction. Default only in
                            direct.
      -e, --excel           For saving results in easy paste to excel format.
                            Default human readable format.