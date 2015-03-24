total_range = 522
segment_size = 20

cluster_match_file_name = "../data/clusters/mm10_clusters_match.txt"
cluster_histogram_file_name = "../data/analysis/clusters_tf_histogram.txt"


def parse_match_file(cluster_file):
    match_list = []
    match_record = {}
    for line in cluster_file:
        if line[0] == '>':
            match_list.append(match_record)
            match_record = {}
        else:
            tf_name = line.split()[0]
            match_record[tf_name] = map(lambda x: abs(int(x)), line.strip().split()[1:])
    match_list.append(match_record)
    match_list.pop(0)
    return match_list


with open(cluster_match_file_name, 'r') as cluster_match_file:
    match_record_list = parse_match_file(cluster_match_file)

histogram_size = (total_range + segment_size) / segment_size
histogram = [0] * histogram_size

for record in match_record_list:
    for value_list in record.values():
        for value in value_list:
            histogram[value / segment_size] += 1

with open(cluster_histogram_file_name, 'w') as cluster_histogram_file:
    cluster_histogram_file.write('\n'.join(map(str, histogram)))