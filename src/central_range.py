import glob
import os

import lib


directory = '../data/converting/'
file_template = 'mm10_*.bed'

bed_file_name_list = []
for root, dirs, files in os.walk(directory):
    for cur_dir in dirs:
        full_cur_dir = os.path.join(root, cur_dir)
        match_template_list = glob.glob1(full_cur_dir, file_template)
        for match_template in match_template_list:
            if '_ucsc' not in match_template:
                full_match_template_path = os.path.join(full_cur_dir, match_template)
                bed_file_name_list.append(full_match_template_path)

for path_bed_file in bed_file_name_list:
    with open(path_bed_file, 'r') as bed_file:
        bed_file_path, bed_file_basename = os.path.split(path_bed_file)
        bed_file_basename_list = os.path.splitext(bed_file_basename)
        bed_converted_file_name = os.path.join(bed_file_path,
                                               bed_file_basename_list[0] + '_ucsc' + bed_file_basename_list[1])
        with open(bed_converted_file_name, 'w') as bed_converted_file:
            for interval_line in bed_file:
                chromosome, start, end = lib.parse_interval_line(interval_line)
                center = (end + start) / 2
                new_start = center - 100
                new_end = center + 100
                bed_converted_file.write(lib.interval_param_to_str(chromosome, new_start, new_end) + '\n')