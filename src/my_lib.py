def searching_result_to_str(interval, results, sequence):
    result_str = ""
    if len(results[0]) == 0:
        return result_str
    else:
        result_str += ">" + interval + '\n'
        for result in results:
            pos_str = ""
            for (pos, score) in result:
                if pos < 0:
                    pos_str += '-'
                    pos += len(sequence)
                pos_str += str(pos) + ' '
            result_str += pos_str + '\n'
    return result_str