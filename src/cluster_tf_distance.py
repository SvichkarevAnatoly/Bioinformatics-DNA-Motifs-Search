cluster_match_file_name = "../data/clusters/mm10_clusters_match.txt"
cluster_tf_dist_file_name = "../data/analysis/clusters_tf_dist.txt"

total_range = 50
segment_size = 1

tf_name_1 = "NANOG"
tf_name_2 = "OCT4"


def parse_match_file(cluster_file):
    match_list = []
    match_record = {}
    for line in cluster_file:
        if line[0] == '>':
            match_list.append(match_record)
            match_record = {}
        else:
            tf_name = line.split()[0]
            match_record[tf_name] = map(int, line.strip().split()[1:])
    match_list.append(match_record)
    match_list.pop(0)
    return match_list


with open(cluster_match_file_name, 'r') as cluster_match_file:
    match_record_list = parse_match_file(cluster_match_file)

tf_1_values_list = []
tf_2_values_list = []
for record in match_record_list:
    tf_1_values_list.append(record[tf_name_1])
    tf_2_values_list.append(record[tf_name_2])

tf_dist_list = []


def tf_min_dist(list1, list2):
    if not list1 or not list2:
        return -1
    else:
        min_dist = abs(list1[0] - list2[0])
        for val1 in list1:
            for val2 in list2:
                d = abs(val1 - val2)
                min_dist = min(d, min_dist)
        return min_dist


for tf_1_values, tf_2_values in zip(tf_1_values_list, tf_2_values_list):
    min_distance = tf_min_dist(tf_1_values, tf_2_values)
    if min_distance >= 0:
        tf_dist_list.append(min_distance)

histogram_size = (total_range + segment_size) / segment_size
histogram = [0] * histogram_size

for dist in tf_dist_list:
    segment_id = dist / segment_size
    if segment_id < histogram_size:
        histogram[segment_id] += 1

with open(cluster_tf_dist_file_name, 'w') as cluster_tf_dist_file:
    cluster_tf_dist_file.write('\n'.join(map(str, histogram)))