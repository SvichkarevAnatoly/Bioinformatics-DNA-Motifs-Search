import re


def pattern_matching_list(text, pattern):
    return [m.start() for m in re.finditer('(?=' + pattern + ')', text)]