def searching_result_to_str(sequence, consensus, results, max_score, seq_number):
    if len(results[0]) == 0:
        return ""
    result_str = str(seq_number) + ":\n"
    result_str += ' ' * 50 + '|' + ' ' * (len(sequence) - 101) + '|' + '\n'
    result_str += sequence + '\n'
    for result in results:
        for (pos, score) in result:
            result_str += ' ' * pos + consensus + \
                          ' ' * (len(sequence) - pos - len(consensus)) + \
                          " Score: " + str(int(score)) + "/" + str(int(max_score)) + '\n'
    return result_str