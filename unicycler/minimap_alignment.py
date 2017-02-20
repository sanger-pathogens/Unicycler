"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functionality related to minimap, which Unicycler uses to seed long read
alignments.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

from .misc import get_nice_header, dim
from collections import defaultdict
from . import log


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
            log.log(dim(line), 3)
            alignment = MinimapAlignment(line, read_dict, ref_dict)
            minimap_alignments[alignment.read_name].append(alignment)
        except (IndexError, ValueError):
            pass
    return minimap_alignments
