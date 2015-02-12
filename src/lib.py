import re


def parse_interval_line(interval_line):
    chr_name, start, end = re.split(":|-", interval_line.strip())
    start = int(start)
    end = int(end)
    return [chr_name, start, end]


def interval_param_to_str(interval_param):
    return interval_param[0] + ':' + str(interval_param[1]) + '-' + str(interval_param[2])


def searching_result_to_str(interval, results, rev_results, sequence_length):
    result_str = ""
    if len(results[0]) == 0 and len(rev_results[0]) == 0:
        return result_str
    else:
        result_str += ">" + interval + '\n'
        for result in results:
            pos_str = ""
            for (pos, score) in result:
                if pos < 0:
                    pos_str += '-'
                    pos += sequence_length
                pos_str += str(pos) + ' '
            result_str += pos_str + '\n'
        for result in rev_results:
            pos_str = ""
            for (pos, score) in result:
                if pos < 0:
                    pos_str += '-'
                    pos += sequence_length
                pos_str += str(pos) + ' '
            result_str += pos_str + '\n'
    return result_str