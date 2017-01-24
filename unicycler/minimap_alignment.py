"""
Class for simple minimap alignments and related functions.

Author: Ryan Wick
email: rrwick@gmail.com
"""

from .misc import get_nice_header
from collections import defaultdict


class MinimapAlignment(object):

    def __init__(self, minimap_line, read_dict, ref_dict):
        line_parts = minimap_line.split('\t')

        self.read_name = line_parts[0]
        self.read_start = int(line_parts[1])
        self.read_end = int(line_parts[2])
        self.read_strand = line_parts[3]

        self.ref_name = get_nice_header(line_parts[4])
        self.ref_start = int(line_parts[5])
        self.ref_end = int(line_parts[6])

        self.read = read_dict[self.read_name]
        self.ref = ref_dict[self.ref_name]

    def get_concise_string(self):
        return ','.join([str(x) for x in [self.read_start, self.read_end, self.read_strand,
                                          self.ref_name, self.ref_start, self.ref_end]])


def line_iterator(string_with_line_breaks):
    """Iterates over a string containing line breaks, one line at a time."""
    prev_newline = -1
    while True:
        next_newline = string_with_line_breaks.find('\n', prev_newline + 1)
        if next_newline < 0:
            break
        yield string_with_line_breaks[prev_newline + 1:next_newline]
        prev_newline = next_newline


def load_minimap_alignments(minimap_alignments_str, read_dict, ref_dict):
    minimap_alignments = defaultdict(list)
    for line in line_iterator(minimap_alignments_str):
        try:
            alignment = MinimapAlignment(line, read_dict, ref_dict)
            minimap_alignments[alignment.read_name].append(alignment)
        except (IndexError, ValueError):
            pass
    return minimap_alignments
