[![Build Status](https://travis-ci.org/SvichkarevAnatoly/Bioinformatics-DNA-Motifs-Search.svg?branch=master)](https://travis-ci.org/SvichkarevAnatoly/Bioinformatics-DNA-Motifs-Search)
[![Code Health](https://landscape.io/github/SvichkarevAnatoly/Bioinformatics-DNA-Motifs-Search/master/landscape.svg?style=flat)](https://landscape.io/github/SvichkarevAnatoly/Bioinformatics-DNA-Motifs-Search/master)

# Bioinformatics-DNA-Motifs-Search
## Requirements
+ [Python 2.7](http://www.python.org)
+ [Biopython](http://biopython.org/)
+ [MOODS](http://www.cs.helsinki.fi/group/pssmfind/)

## Contain utilities
+ [bed_center_extender.py](#usage-bed_center_extenderpy--h--l-length--o-outfile-bedfile)
+ [pattern_matching.py](#usage-pattern_matchingpy--h--o-output--tf-tf-tf---th-threshold--rc--b--e-fasta-pwm)
+ [excel_ucsc_ids_to_bed_order.py](#usage-excel_ucsc_ids_to_bed_orderpy--h--b-bed--o-output-excel)
+ [pwm_generator.py](#usage-pwm_generatorpy--h--o-output-seqs)
+ [logo_generator.py](#usage-logo_generatorpy--h--s-seqs---p-pwm--o-output)

## Usage
#### usage: bed_center_extender.py \[-h] \[-l LENGTH] \[-o OUTFILE] bedfile

    Central extension each interval of the specified file to the same length

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

#### usage: pattern_matching.py \[-h] \[-o \[OUTPUT]] \[-tf TF \[TF ...]] \[-th THRESHOLD] \[-rc] \[-b] \[-e] fasta pwm
    
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
                            reported at position [position - sequence_length] in
                            complement of input sequence, which is always
                            negative. The actual hit site for any hit is always
                            seq[pos, pos + matrix_length]. Default False.
      -e, --excel           For saving results in easy paste to excel format.
                            Default human readable format.
                            
#### usage: excel_ucsc_ids_to_bed_order.py \[-h] \[-b \[BED]] \[-o \[OUTPUT]] excel
    
    Converting matching results in plain text excel format in bed file orders
    
    positional arguments:
      excel                 text file with excel matching
    
    optional arguments:
      -h, --help            show this help message and exit
      -b [BED], --bed [BED]
                            bed file with intervals. Order identifiers to order in
                            bed file. Default not order.
      -o [OUTPUT], --output [OUTPUT]
                            output file with formatted matching results. Default
                            stdout.
                            
#### usage: pwm_generator.py \[-h] \[-o \[OUTPUT]] seqs

    Create position weight matrices (PWM) from DNA sequences
    
    positional arguments:
      seqs                  file with DNA sequences same length. One sequence in
                            on line
    
    optional arguments:
      -h, --help            show this help message and exit
      -o [OUTPUT], --output [OUTPUT]
                            output file with PWM

#### usage: logo_generator.py \[-h] (-s SEQS | -p PWM) -o OUTPUT

    Create sequence logo from sequences
    
    optional arguments:
      -h, --help            show this help message and exit
      -s SEQS, --seqs SEQS  file with DNA sequences same length. One sequence in
                            on line
      -p PWM, --pwm PWM     file with PWM.
      -o OUTPUT, --output OUTPUT
                            output file with logo in vector SVG format
