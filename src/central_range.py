import glob
import os

directory = '../data/converting/'
file_template = 'mm10_*.bed'

bed_file_name_list = []
for root, dirs, files in os.walk(directory):
    for cur_dir in dirs:
        full_cur_dir = os.path.join(root, cur_dir)
        match_template_list = glob.glob1(full_cur_dir, file_template)
        for match_template in match_template_list:
            full_match_template_path = os.path.join(full_cur_dir, match_template)
            bed_file_name_list.append(full_match_template_path)

for path_bed_file in bed_file_name_list:
    print path_bed_file