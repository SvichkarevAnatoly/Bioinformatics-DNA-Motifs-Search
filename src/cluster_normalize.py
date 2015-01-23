import os
import re

interval_range = 0
bed_file_name = "../data/clusters/mm10_clusters.bed"

with open(bed_file_name, 'r') as bed_file_tmp:
    for interval_line in bed_file_tmp:
        chromosome, start, end = re.split(":|-", interval_line.strip())
        start = int(start)
        end = int(end)
        diff = end - start
        interval_range = max(interval_range, diff)
interval_range += 1

with open(bed_file_name, 'r') as bed_file:
    bed_file_path, bed_file_basename = os.path.split(bed_file_name)
    bed_file_basename_list = os.path.splitext(bed_file_basename)
    bed_converted_file_name = os.path.join(bed_file_path,
                                           bed_file_basename_list[0] + '_ucsc' + bed_file_basename_list[1])
    with open(bed_converted_file_name, 'w') as bed_converted_file:
        for interval_line in bed_file:
            chromosome, start, end = re.split(":|-", interval_line.strip())
            start = int(start)
            end = int(end)
            center = (end + start) / 2
            new_start = str(center - interval_range / 2)
            new_end = str(center + interval_range / 2)
            new_interval_line = chromosome + ' ' + new_start + ' ' + new_end
            bed_converted_file.write(new_interval_line + '\n')