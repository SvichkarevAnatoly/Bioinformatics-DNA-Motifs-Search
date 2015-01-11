#!/bin/bash

input_file="$1"
extension_input_file="${input_file##*.}"
name_input_file="${input_file%.*}"

if [ "$#" = 1 ]; then
    output_file="$name_input_file"_out."$extension_input_file"
    sed 's/[:-]/ /g' "$input_file" > "$output_file"
    echo Converted: "$input_file"
else
    output_file=$2
    sed 's/[:-]/ /g' "$input_file" > "$output_file"
    echo Converted: "$input_file" to "$output_file"
fi