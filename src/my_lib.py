import MOODS


def pattern_matching_list(text, pattern):
    # TATA
    matrix = [[0, 1, 0, 1],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [1, 0, 1, 0]]

    results = MOODS.search(text, [matrix], 3, convert_log_odds=False,
                           threshold_from_p=False)
    return list(results)


def matching_list_to_string(matching_list):
    string = ''
    for i in matching_list:
        for (position, score) in i:
            string += "Position: " + str(position) + " Score: " + str(score) + '\n'
    return string