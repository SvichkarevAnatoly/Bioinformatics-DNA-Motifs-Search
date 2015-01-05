import MOODS


def pattern_to_pwm(pattern):
    """
    Write a pattern sequence as a position weight matrix.
    """
    # Create 4 lists of length equal to primer's length.
    matrix = [[0] * len(pattern) for i in range(4)]
    # List of correspondance IUPAC.
    IUPAC = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "U": ["U"],
        "R": ["G", "A"],
        "Y": ["T", "C"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"]
    }
    # Position of nucleotides in the PWM.
    dico = {"A": 0, "C": 1, "G": 2, "T": 3}
    # Read each IUPAC letter in the primer.
    for index, letter in enumerate(pattern):
        for nuc in IUPAC.get(letter):
            i = dico.get(nuc)
            matrix[i][index] = 1
    return matrix


def compute_thresholds(pattern):
    thresholds = 0.8 * len(pattern)
    if len(pattern) - thresholds < 1:
        thresholds = len(pattern) - 1
    return thresholds


def pattern_matching_list(text, pattern):
    thresholds = compute_thresholds(pattern)
    matrix = pattern_to_pwm(pattern)
    results = MOODS.search(text, [matrix], thresholds, convert_log_odds=False,
                           threshold_from_p=False)
    return list(results)


def matching_list_to_string(matching_list):
    string = ''
    for i in matching_list:
        for (position, score) in i:
            string += "Position: " + str(position) + " Score: " + str(score) + '\n'
    return string