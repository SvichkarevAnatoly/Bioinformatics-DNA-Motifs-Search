#!/bin/bash

echo bed center extender
python ../../../../src/utils/bed_center_extender.py mm10.bed -l 500 -o mm10_500.bed

echo bed to UCSC Table Browser
./../../../../scripts/bed_to_UCSC_Table_Browser.sh mm10_500.bed mm10_500_ucsc.bed
