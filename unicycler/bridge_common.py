"""
Bridges are links between two single-copy segments in an assembly graph.
Bridges can come from multiple sources, so there are a few different classes which all share
the important members/methods (for duck-typing purposes).

This module has a few functions that are common to multiple bridge types.

Author: Ryan Wick
email: rrwick@gmail.com
"""


import math
from .misc import weighted_average


def get_mean_depth(seg_1, seg_2, graph):
    """
    Returns the mean depth of the two segments, weighted by their length.
    """
    return weighted_average(seg_1.depth, seg_2.depth,
                            seg_1.get_length_no_overlap(graph.overlap),
                            seg_2.get_length_no_overlap(graph.overlap))


def get_bridge_str(bridge):
    """
    Returns a bridge sequence in human-readable form.
    """
    bridge_str = str(bridge.start_segment) + ' -> '
    if bridge.graph_path:
        bridge_str += ', '.join([str(x) for x in bridge.graph_path]) + ' -> '
    bridge_str += str(bridge.end_segment)
    return bridge_str


def get_depth_agreement_factor(start_seg_depth, end_seg_depth):
    """
    This function is set up such that:
      * equal depths return 1.0
      * similar depths return a value near 1.0
      * more divergent depths return a much lower value:
          a ratio of 1.35 return a value of about 0.5
          a ratio of 2.06 return a value of about 0.1
      * very different depths return a value near 0.0
    https://www.desmos.com/calculator
        y=\frac{1}{1+10^{2\left(\log \left(x-1\right)+0.45\right)}}
        y=\frac{1}{1+10^{2\left(\log \left(\frac{1}{x}-1\right)+0.45\right)}}
    """
    larger_depth = max(start_seg_depth, end_seg_depth)
    smaller_depth = min(start_seg_depth, end_seg_depth)
    if larger_depth == 0.0 or smaller_depth == 0.0:
        return 0.0
    elif larger_depth == smaller_depth:
        return 1.0
    else:
        ratio = larger_depth / smaller_depth
        return 1.0 / (1.0 + 10.0 ** (2 * (math.log10(ratio - 1.0) + 0.45)))
