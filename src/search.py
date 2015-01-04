import sys

from my_lib import pattern_matching_list

input_file_str = sys.argv[1]
with open(input_file_str, 'r') as input_file:
    genome = input_file.readline().strip()
    pattern = input_file.readline().strip()

matching_list = pattern_matching_list(genome, pattern)

if len(sys.argv) == 3:
    with open(sys.argv[2], 'w') as output_file:
        output_file.write(' '.join(map(str, matching_list)) + '\n')
else:
    print ' '.join(map(str, matching_list))
