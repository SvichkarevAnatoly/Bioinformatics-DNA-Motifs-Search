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
    cluster_match_list = parse_match_file(cluster_match_file)

pass