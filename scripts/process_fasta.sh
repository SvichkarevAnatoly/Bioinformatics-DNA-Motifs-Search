#!/bin/bash

tf=$(basename $PWD)
echo TF: ${tf}

echo pattern mathing
python ../../../../src/utils/pattern_matching.py mm10_500.fa ../13TF_pwm.dat -tf ${tf} -rc -e -o mm10_500_excel.txt

echo best matches
awk '{print $NF}' mm10_500_excel.txt | sed "/repeatMasking=none]/d" > mm10_500_best_matches.txt

echo pwm generator
python ../../../../src/utils/pwm_generator.py mm10_500_best_matches.txt -m ${tf} -o pwm.dat

echo logo generator
python ../../../../src/utils/logo_generator.py -p pwm.dat -o ${tf}.svg
