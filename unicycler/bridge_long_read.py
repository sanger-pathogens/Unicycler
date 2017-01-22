"""
Long read bridges are the big important type of bridge for a hybrid Unicycler assembly. They are
made using long reads which align to multiple segments in the graph.

Author: Ryan Wick
email: rrwick@gmail.com
"""

from multiprocessing.dummy import Pool as ThreadPool
import time
import math
import statistics
from collections import defaultdict
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor
from .misc import float_to_str, reverse_complement, flip_number_order, score_function, \
    print_table, get_right_arrow
from .cpp_function_wrappers import multiple_sequence_alignment
from . import settings
from .path_finding import get_best_paths_for_seq


class LongReadBridge(object):
    """
    This class describes a bridge created from long read alignments.
    """
    def __init__(self, graph, start, end):

        # The numbers of the two single-copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The individual reads contributing to the bridge. The sequences/qualities are not for the
        # entire read, just the part in the bridge. In the case of overlapping alignments, this has
        # the overlap size instead of the sequence.
        self.full_span_reads = []

        # The bridge can also contain incomplete read sequences which don't bridge the entire span
        # between the start and end segments. These are still useful as they can contribute to the
        # consensus sequence.
        self.start_only_reads = []
        self.end_only_reads = []

        # The consensus of all read sequences. If there is only one read, then this is the same as
        # that read's sequence. If there are multiple reads, this is hopefully a more accurate
        # sequence. If this is a case where the start and end segments overlap, then there will not
        # be a consensus sequence and consensus_length will be the mean of read_lengths (a negative
        # value).
        self.consensus_sequence = ''

        # The path through the unbridged graph, if one was found.
        self.graph_path = []
        self.all_paths = []

        # The bridge sequence, gotten from the graph path if a good path was found. Otherwise it's
        # from the consensus sequence.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        # A score used to determine the order of bridge application.
        self.quality = 1.0

        # Whether or not the long read bridge got its sequence from a graph path (True) or from the
        # read consensus (False).
        self.path_support = None

        # When a bridge is applied, the segments in the bridge may have their depth reduced
        # accordingly. This member stores which segments have had their depth reduced and by how
        # much due to this bridge's application. It is stored so if this bridge is later deleted,
        # we can restore the depth to the segments.
        self.segments_reduced_depth = []

        self.graph = graph

    def __repr__(self):
        return 'long read bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    def predicted_time_to_finalise(self):
        """
        This function very roughly predicts how long the bridge will take to finalise. It's not
        meant to be particularly accurate, but can hopefully be used to roughly order the bridges
        from slow to fast.
        """
        total_full_span_seq_length = 0
        full_span_seq_count = 0
        for full_span in self.full_span_reads:
            if not isinstance(full_span[0], int):
                total_full_span_seq_length += len(full_span[0])
                full_span_seq_count += 1
        if not full_span_seq_count:
            mean_full_span_seq_length = 0.0
        else:
            mean_full_span_seq_length = total_full_span_seq_length / full_span_seq_count

        if full_span_seq_count > 1:
            predicted_consensus_time = ((1.34e-9 * (total_full_span_seq_length ** 2)) +
                                        (2.76e-5 * total_full_span_seq_length))
        else:
            predicted_consensus_time = 0.0
        predicted_path_time = ((1.78e-7 * (mean_full_span_seq_length ** 2)) +
                               (3.75e-3 * mean_full_span_seq_length))

        return predicted_consensus_time + predicted_path_time

    def finalise(self, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size,
                 expected_linear_seqs):
        """
        Determines the consensus sequence for the bridge, attempts to find it in the graph and
        assigns a quality score to the bridge. This is the big performance-intensive step of long
        read bridging!
        """
        start_seg = self.graph.segments[abs(self.start_segment)]
        end_seg = self.graph.segments[abs(self.end_segment)]

        output = [str(self.start_segment), str(self.end_segment), str(len(self.full_span_reads))]

        start_alignment_scaled_scores = [x[2].scaled_score for x in self.full_span_reads]
        end_alignment_scaled_scores = [x[3].scaled_score for x in self.full_span_reads]
        best_overall_scaled_score = min(max(start_alignment_scaled_scores),
                                        max(end_alignment_scaled_scores))
        alignment_scaled_scores = start_alignment_scaled_scores + end_alignment_scaled_scores
        mean_alignment_scaled_score = statistics.mean(alignment_scaled_scores)
        read_to_ref_ratios = [x[2].get_read_to_ref_ratio() for x in self.full_span_reads] + \
                             [x[3].get_read_to_ref_ratio() for x in self.full_span_reads]
        mean_read_to_ref_ratio = statistics.mean(read_to_ref_ratios)

        # Partition the full span reads into two groups: those with negative numbers (implying that
        # the two segments overlap) and those with actual sequences.
        full_spans_without_seq = []
        full_spans_with_seq = []
        for full_span in self.full_span_reads:
            if isinstance(full_span[0], int):
                full_spans_without_seq.append(full_span)
            else:
                full_spans_with_seq.append(full_span)

        # There shouldn't usually be both full spans with sequence and without. If there are some
        # of each, we'll throw out the minority group.
        if full_spans_with_seq and full_spans_without_seq:
            if len(full_spans_without_seq) > len(full_spans_with_seq):
                full_spans_with_seq = []
            else:
                full_spans_without_seq = []

        # For full spans with sequence, we perform a MSA and get a consensus sequence.
        if full_spans_with_seq:

            # Full-span sequences are faster in the MSA and are very useful for the consensus
            # sequence, so we want to use a lot! But we still set an upper limit for cases of very
            # high read depth, as getting the consensus sequence will be slow if there are too many
            # reads.
            max_full_span = settings.MAX_READS_FOR_CONSENSUS
            if len(full_spans_with_seq) <= max_full_span:
                full_span_seqs = [x[0] for x in full_spans_with_seq]
                full_span_quals = [x[1] for x in full_spans_with_seq]
            else:
                sorted_full_span = sorted(full_spans_with_seq,
                                          key=lambda x: x[2].raw_score + x[3].raw_score,
                                          reverse=True)
                full_span_seqs = [x[0] for x in sorted_full_span[:max_full_span]]
                full_span_quals = [x[1] for x in sorted_full_span[:max_full_span]]

            # # Start-only and end-only sequences make the MSA a lot slower, contribute less to the
            # # consensus and can screw up the MSA if they are too abundant. So if there are too
            # # many of them, we only include the ones with the best alignments.
            # # Also, if we have a lot of full span reads, there's no need to include any partial
            # # sequences at all.
            # max_partial = min(3, len(full_span_seqs))  # TO DO: make this a parameter?
            # if len(full_span_seqs) > max_full_span / 2:
            #     max_partial = 0

            # TEMP - turning partial sequences off to test the effect. If I decided to turn them
            # permanently off, I can simplify quite a lot of code.
            max_partial = 0

            if len(self.start_only_reads) <= max_partial:
                start_only_seqs = [x[0] for x in self.start_only_reads]
                start_only_quals = [x[1] for x in self.start_only_reads]
            else:
                sorted_start_only = sorted(self.start_only_reads,
                                           key=lambda x: x[2].raw_score, reverse=True)
                start_only_seqs = [x[0] for x in sorted_start_only[:max_partial]]
                start_only_quals = [x[1] for x in sorted_start_only[:max_partial]]

            if len(self.end_only_reads) <= max_partial:
                end_only_seqs = [x[0] for x in self.end_only_reads]
                end_only_quals = [x[1] for x in self.end_only_reads]
            else:
                sorted_end_only = sorted(self.end_only_reads,
                                         key=lambda x: x[2].raw_score, reverse=True)
                end_only_seqs = [x[0] for x in sorted_end_only[:max_partial]]
                end_only_quals = [x[1] for x in sorted_end_only[:max_partial]]

            consensus_start_time = time.time()
            self.consensus_sequence = multiple_sequence_alignment(full_span_seqs, full_span_quals,
                                                                  start_only_seqs, start_only_quals,
                                                                  end_only_seqs, end_only_quals,
                                                                  scoring_scheme)[0]
            consensus_time = time.time() - consensus_start_time

            output.append(str(len(self.consensus_sequence)))
            output.append(float_to_str(consensus_time, 1))
            # if verbosity > 3:
            #     output += '                           ' + self.consensus_sequence + '\n'

            # We now make an expected scaled score for an alignment between the consensus and a
            # graph path. I.e. when we find a path in the graph for this consensus, this is about
            # how well it should align. This goes up with alignment scores (better alignment scores
            # imply better read sequences) and read count (more reads can give a better consensus).
            #     https://www.desmos.com/calculator
            #     y=100\cdot \left(\left(1-\frac{a}{100}\right)
            #       \left(1-\frac{3}{2+x}\right)+\frac{a}{100}\right)
            #     y = expected_scaled_score, x = num_span_reads, a = mean_alignment_scaled_score
            num_span_reads = len(self.full_span_reads)
            expected_scaled_score = 100.0 * ((1.0 - mean_alignment_scaled_score / 100.0) *
                                             (1.0 - (3.0 / (2.0 + num_span_reads))) +
                                             mean_alignment_scaled_score / 100.0)
            expected_scaled_score = max(expected_scaled_score, best_overall_scaled_score)

            # We can also predict the ratio in length between the consensus sequence and the
            # graph path. For a low number of input reads, this is close to the mean ratio for
            # the read alignments, but it approaches 1 as we have more reads and expect our
            # consensus to be better.
            #     https://www.desmos.com/calculator
            #     y=\left(a-1\right)\left(\frac{b}{x+b-1}\right)+1
            expected_consensus_to_ref_ratio = 1.0 + (mean_read_to_ref_ratio - 1.0) * \
                                                    (4 / (4 + num_span_reads - 1))
            target_path_length = (2 * self.graph.overlap) + \
                int(round((len(self.consensus_sequence) / expected_consensus_to_ref_ratio)))

            mean_overlap = 0

        # For full spans without sequence, we simply need a mean distance.
        else:
            self.consensus_sequence = ''
            mean_overlap = int(round(sum(x[0] for x in full_spans_without_seq) /
                                     len(full_spans_without_seq)))
            output.append(str(mean_overlap))
            output.append('')
            target_path_length = mean_overlap + (2 * self.graph.overlap)
            expected_scaled_score = 100.0

        output.append(str(target_path_length))

        path_start_time = time.time()
        self.all_paths, progressive_path_search = \
            get_best_paths_for_seq(self.graph, self.start_segment, self.end_segment,
                                   target_path_length, self.consensus_sequence, scoring_scheme,
                                   expected_scaled_score)
        path_time = time.time() - path_start_time

        output.append(str(len(self.all_paths)))
        output.append('progressive' if progressive_path_search else 'exhaustive')
        output.append(float_to_str(path_time, 1))

        # If paths were found, use a path sequence for the bridge.
        if self.all_paths:
            if self.all_paths[0][0]:
                output.append(', '.join(str(x) for x in self.all_paths[0][0]))
            else:
                output.append('direct connection')
            best_path_len = self.graph.get_bridge_path_length(self.all_paths[0][0])
            output.append(str(best_path_len))
            output.append(float_to_str(self.all_paths[0][1], 1))
            output.append(float_to_str(self.all_paths[0][3], 2))
            output.append(str(self.all_paths[0][2]))

            self.graph_path = self.all_paths[0][0]
            self.bridge_sequence = self.graph.get_bridge_path_sequence(self.graph_path,
                                                                       self.start_segment)
            self.path_support = True

            # We start this bridge's quality using a function that takes into account the
            # actual, expected and minimum acceptable scores. If the actual scaled score is 100,
            # this function gives 1. If it is the expected value, this function gives a number
            # around 0.7. Then as the actual score approaches the minimum acceptable score, this
            # function approaches 0.
            #     https://www.desmos.com/calculator
            #     y=\left(\frac{1}{1+2^{a-x}}\right)^{0.5}
            #     y = self.quality, x = actual_scaled_score, a = expected_scaled_score,
            #     b = min_acceptable_scaled_score
            actual_scaled_score = self.all_paths[0][3]
            self.quality = math.sqrt(1.0 /
                                     (1.0 + 2.0 ** (expected_scaled_score - actual_scaled_score)))

        # If a path wasn't found, the consensus sequence is the bridge (with the overlaps added).
        else:
            self.graph_path = []
            self.path_support = False
            output += ['', '', '', '', '']
            start_overlap = \
                self.graph.seq_from_signed_seg_num(self.start_segment)[-self.graph.overlap:]
            end_overlap = \
                self.graph.seq_from_signed_seg_num(self.end_segment)[:self.graph.overlap]

            # If there is a consensus sequence, we simply tack the overlaps onto its ends.
            if self.consensus_sequence:
                self.bridge_sequence = start_overlap + self.consensus_sequence + end_overlap

            # If there is the consensus sequence is exactly zero, then the two segments butt up
            # exactly.
            elif mean_overlap == 0:
                self.bridge_sequence = start_overlap + end_overlap

            # If the consensus sequence is negative, that implies the start and the end segments
            # overlap. In this case we need to find the exact overlap. If we don't find it, we
            # just glue the sequences together, like we did for mean_overlap == 0.
            else:
                larger_overlaps = list(range(abs(mean_overlap), self.graph.overlap))
                smaller_overlaps = list(range(abs(mean_overlap) - 1, 0, -1))
                test_overlaps = []
                for i in range(max(len(larger_overlaps), len(smaller_overlaps))):
                    if i < len(larger_overlaps):
                        test_overlaps.append(larger_overlaps[i])
                    if i < len(smaller_overlaps):
                        test_overlaps.append(smaller_overlaps[i])
                for test_overlap in test_overlaps:
                    if start_overlap[-test_overlap:] == end_overlap[:test_overlap]:
                        actual_overlap = test_overlap
                        break
                else:
                    actual_overlap = 0
                self.bridge_sequence = start_overlap + end_overlap[actual_overlap:]

            # The quality of non-graph-path-supported bridges depends on the number of dead ends
            # and whether or not linear sequences are expected.
            dead_end_count = 0
            if self.graph.ends_with_dead_end(self.start_segment):
                dead_end_count += 1
            if self.graph.starts_with_dead_end(self.end_segment):
                dead_end_count += 1

            if expected_linear_seqs:
                if dead_end_count == 2:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS_WITH_LINEAR_SEQS
                elif dead_end_count == 1:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_ONE_DEAD_END_WITH_LINEAR_SEQS
                else:  # dead_end_count == 0
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS_WITH_LINEAR_SEQS
            else:
                if dead_end_count == 2:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS
                elif dead_end_count == 1:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_ONE_DEAD_END
                else:  # dead_end_count == 0
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS

        # Expected read count is determined using the read lengths and bridge size. For a given
        # read length and bridge, there are an estimable number of positions where a read of that
        # length would be able to contribute to the bridge. This is used to get the probability
        # that any read would create a bridge, and totalling those up gives us our estimated count.
        bridge_len = max(len(self.bridge_sequence), self.graph.overlap)
        min_read_len = (2 * min_alignment_length) + bridge_len - (2 * self.graph.overlap)
        total_possible_placements = 0
        for read_len, count in read_lengths.items():
            if read_len < min_read_len:
                continue
            possible_read_placements = read_len - min_read_len + 1
            possible_read_placements *= count
            possible_read_placements *= max(self.depth, 1)
            total_possible_placements += possible_read_placements
        expected_read_count = total_possible_placements / estimated_genome_size
        actual_read_count = len(self.full_span_reads)

        # Adjust the expected read count down, especially for higher values.
        # TO DO: reevaluate this step - is it necessary?
        expected_read_count = reduce_expected_count(expected_read_count, 30, 0.5)

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # as it implies that they aren't actually single-copy or on the same piece of DNA.
        depth_agreement_factor = get_depth_agreement_factor(start_seg.depth, end_seg.depth)
        self.quality *= depth_agreement_factor

        # The number of reads which contribute to a bridge is a big deal, so the read count factor
        # scales linearly. This value is capped at 1, which means that bridges with too few reads
        # are punished but bridges with excess reads are not rewarded.
        read_count_factor = min(1.0, actual_read_count / expected_read_count)
        self.quality *= read_count_factor

        # The length of alignments to the start/end segments is positively correlated with quality
        # to reward bridges with long alignments. Specifically, we want there to be at least one
        # spanning read with a long alignment to the start segment and at least one spanning read
        # with a long alignment to the end segment.
        longest_start_alignment = max(x[2].get_aligned_ref_length() for x in self.full_span_reads)
        longest_end_alignment = max(x[3].get_aligned_ref_length() for x in self.full_span_reads)
        alignment_length = min(longest_start_alignment, longest_end_alignment)
        align_length_factor = score_function(alignment_length, min_alignment_length * 4)
        self.quality *= align_length_factor

        # The mean alignment score to the start/end segments is positively correlated with quality,
        # so bridges with high quality alignments are rewarded. Specifically, we want there to be at
        # least one spanning read with a high quality alignment to the start segment and at least
        # one spanning read with a high quality alignment to the end segment.
        best_start_alignment = max(x[2].scaled_score for x in self.full_span_reads)
        best_end_alignment = max(x[3].scaled_score for x in self.full_span_reads)
        alignment_quality = min(best_start_alignment, best_end_alignment)
        align_score_factor = alignment_quality / 100.0
        self.quality *= align_score_factor

        # Bridges between long start/end segments are rewarded, as they are more likely to actually
        # be single-copy. We apply a length factor for both the start and the end segments,
        # and then apply the smaller of two again. This is to punish cases where both segments
        # are not long.
        start_length_factor = score_function(start_seg.get_length(), min_alignment_length * 4)
        self.quality *= start_length_factor
        end_length_factor = score_function(end_seg.get_length(), min_alignment_length * 4)
        self.quality *= end_length_factor
        smaller_length_factor = min(start_length_factor, end_length_factor)
        self.quality *= smaller_length_factor

        # We finalise the quality to a range of 0 to 100. We also use the sqrt function to pull
        # the scores up a bit (otherwise they tend to hang near the bottom of the range).
        self.quality = 100.0 * math.sqrt(self.quality)

        output.append(self.quality)

        return output

    def set_path_based_on_availability(self, graph, unbridged_graph):
        """
        This function will change a bridge's graph path based on what's currently available. This
        is to handle the case where a bridge has multiple possible graph paths, but its first
        choice isn't available anymore (because other bridges used up those segments).
        It has to balance the path quality with the path availability to make a choice.
        """
        best_path = self.all_paths[0][0]
        best_sequence = unbridged_graph.get_bridge_path_sequence(best_path, self.start_segment)
        best_scaled_score = self.all_paths[0][3]
        best_availability = graph.get_path_availability(best_path)
        for i in range(1, len(self.all_paths)):
            potential_path = self.all_paths[i][0]
            potential_scaled_score = self.all_paths[i][3]
            potential_availability = graph.get_path_availability(potential_path)

            # relative_score measures how much worse this path aligned than the current best.
            # Differences matter more close to scores of 100. E.g. 99 to 97 is a big drop in
            # score, but 79 to 77 is not as large of a drop.
            if potential_scaled_score == 100.0:
                relative_score = 1.0
            else:
                relative_score = (100.0 - best_scaled_score) / (100.0 - potential_scaled_score)
                relative_score = min(1.0, relative_score)

            # relative_availability measures how much more available this path is than the current
            # best. We use 1.1 (instead of 1.0) in the equation to attenuate the affect of
            # availability a bit (because score is more important).
            relative_availability = (1.1 - best_availability) / (1.1 - potential_availability)
            relative_availability = min(2.0, relative_availability)

            # If this path looks better than our current best (considering both score and
            # availability), then it becomes the new best.
            if relative_score * relative_availability > 1.0:
                best_path = potential_path
                best_sequence = unbridged_graph.get_bridge_path_sequence(potential_path,
                                                                         self.start_segment)
                best_scaled_score = potential_scaled_score
                best_availability = potential_availability

        self.graph_path = best_path
        self.bridge_sequence = best_sequence

    def contains_full_span_sequence(self):
        """
        Some LongReadBridge objects bridge two close segments, and therefore do not actually have
        any bridging read sequence (just a negative number implying overlap). This function returns
        True if any full span read sequences exist and False if they are just negative numbers.
        """
        for full_span_read in self.full_span_reads:
            seq_or_num = full_span_read[0]
            if not isinstance(seq_or_num, int):
                return True
        return False

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 2

    @staticmethod
    def get_type_name():
        """
        Returns the of the bridge types.
        """
        return 'long read'


def create_long_read_bridges(graph, read_dict, read_names, single_copy_segments, verbosity,
                             existing_bridges, min_scaled_score, threads, scoring_scheme,
                             min_alignment_length, expected_linear_seqs, min_bridge_qual):
    """
    Makes bridges between single-copy segments using the alignments in the long reads.
    """
    single_copy_seg_num_set = set()
    for seg in single_copy_segments:
        single_copy_seg_num_set.add(seg.number)

    # This dictionary will collect the read sequences which span between two single-copy segments.
    # These are the most useful sequences and will be used to either create a new bridge or enhance
    # an existing bridge.
    # Key = tuple of signed segment numbers (the segments being bridged)
    # Value = list of tuples containing the bridging sequence and the single-copy segment
    #         alignments.
    spanning_read_seqs = defaultdict(list)

    # This dictionary will collect all of the read sequences which don't manage to span between
    # two single-copy segments, but they do overlap one single-copy segment and can therefore be
    # useful for generating consensus sequences in a bridge.
    # Key = signed segment number, Value = list of tuples containing the sequence and alignment
    overlapping_read_seqs = defaultdict(list)

    for read_name in read_names:
        read = read_dict[read_name]
        alignments = get_single_copy_alignments(read, single_copy_seg_num_set, min_scaled_score)
        if not alignments:
            continue

        # If the code got here, then we have some alignments to single-copy segments. We grab
        # neighbouring pairs of alignments, starting with the highest scoring ones and work our
        # way down. This means that we should have a pair for each neighbouring alignment, but
        # potentially also more distant pairs if the alignments are strong.
        already_added = set()
        sorted_alignments = sorted(alignments, key=lambda x: x.raw_score, reverse=True)
        available_alignments = []
        for alignment in sorted_alignments:

            # If the alignment being added is to a reference that has already been added but in the
            # opposite direction, then we don't include it. E.g. we don't add an alignment for 10
            # if we already have an alignment for -10. This is because there's no legitimate way
            # for a single-copy segment to appear in the same read in two different directions. The
            # same direction is okay, as that can happen with a circular piece of DNA, but opposite
            # directions implies multi-copy.
            opposite_num = -alignment.get_signed_ref_num()
            if opposite_num in set(x.get_signed_ref_num() for x in available_alignments):
                continue

            available_alignments.append(alignment)
            available_alignments = sorted(available_alignments,
                                          key=lambda x: x.read_start_positive_strand())
            if len(available_alignments) < 2:
                continue

            for i in range(len(available_alignments)):
                if i < len(available_alignments) - 1:
                    alignment_1 = available_alignments[i]
                    alignment_2 = available_alignments[i + 1]

                # Special case: when the first and last alignments are to the same graph segment,
                # make a bridge for them, even if they aren't a particularly high scoring pair of
                # alignments. This can help to circularise plasmids which are very tied up with
                # other, similar plasmids.
                elif available_alignments[0].ref.name == available_alignments[-1].ref.name:
                    alignment_1 = available_alignments[0]
                    alignment_2 = available_alignments[-1]
                else:
                    continue

                # Standardise the order so we don't end up with both directions (e.g. 5 to -6 and
                # 6 to -5) in spanning_read_seqs.
                seg_nums, flipped = flip_number_order(alignment_1.get_signed_ref_num(),
                                                      alignment_2.get_signed_ref_num())
                if seg_nums not in already_added:
                    bridge_start = alignment_1.read_end_positive_strand()
                    bridge_end = alignment_2.read_start_positive_strand()

                    if bridge_end > bridge_start:
                        bridge_seq = read.sequence[bridge_start:bridge_end]
                        bridge_qual = read.qualities[bridge_start:bridge_end]
                        if flipped:
                            bridge_seq = reverse_complement(bridge_seq)
                            bridge_qual = bridge_qual[::-1]
                    else:
                        bridge_seq = bridge_end - bridge_start  # 0 or a negative number
                        bridge_qual = ''

                    spanning_read_seqs[seg_nums].append((bridge_seq, bridge_qual, alignment_1,
                                                         alignment_2))
                    already_added.add(seg_nums)

        # At this point all of the alignments have been added and we are interested in the first
        # and last alignments (which may be the same if there's only one). If the read extends
        # past these alignments, then the overlapping part might be useful for consensus sequences.
        first_alignment = available_alignments[0]
        start_overlap, start_qual = first_alignment.get_start_overlapping_read_seq()
        if start_overlap:
            seg_num = -first_alignment.get_signed_ref_num()
            seq = reverse_complement(start_overlap)
            qual = start_qual[::-1]
            if seg_num not in overlapping_read_seqs:
                overlapping_read_seqs[seg_num] = []
            overlapping_read_seqs[seg_num].append((seq, qual, first_alignment))
        last_alignment = available_alignments[-1]
        end_overlap, end_qual = last_alignment.get_end_overlapping_read_seq()
        if end_overlap:
            seg_num = last_alignment.get_signed_ref_num()
            overlapping_read_seqs[seg_num].append((end_overlap, end_qual, last_alignment))

    # If a bridge already exists for a spanning sequence, we add the sequence to the bridge. If
    # not, we create a new bridge and add it.
    new_bridges = []
    for seg_nums, span in spanning_read_seqs.items():
        start, end = seg_nums

        # If the start and end are the same and already exclusively connect (i.e. if this segment
        # is already circular), then skip - there's no need to bridge.
        if start == end and graph.get_downstream_seg_nums(start) == [start] and \
                graph.get_upstream_seg_nums(start) == [start]:
            continue

        for existing_bridge in existing_bridges:
            if isinstance(existing_bridge, LongReadBridge) and \
                    existing_bridge.start_segment == start and existing_bridge.end_segment == end:
                matching_bridge = existing_bridge
                break
        else:
            new_bridge = LongReadBridge(graph, start, end)
            new_bridges.append(new_bridge)
            matching_bridge = new_bridge
        matching_bridge.full_span_reads += span
    all_bridges = existing_bridges + new_bridges
    all_bridges = sorted(all_bridges, key=lambda x: (x.start_segment, x.end_segment))

    # Add overlapping sequences to appropriate bridges, but only if they have some full span
    # sequence (if their full span reads only have overlap-indicating negative numbers, then we
    # won't be doing a consensus sequence and don't need overlapping sequences).
    for seg_num, overlaps in overlapping_read_seqs.items():
        for bridge in all_bridges:
            if not isinstance(bridge, LongReadBridge) or not bridge.contains_full_span_sequence():
                continue
            start_overlap = (bridge.start_segment == seg_num)
            end_overlap = (bridge.end_segment == -seg_num)
            if start_overlap or end_overlap:
                for overlap in overlaps:
                    if start_overlap:
                        bridge.start_only_reads.append(overlap)
                    elif end_overlap:
                        overlap_seq, overlap_qual, alignment = overlap
                        bridge.start_only_reads.append((reverse_complement(overlap_seq),
                                                        overlap_qual[::-1], alignment))

    # During finalisation, we will compare the expected read count to the actual read count for
    # each bridge. To do this, we'll need the lengths of all reads (excluding those with no
    # alignments). We also need an estimate of the genome size.
    read_lengths = defaultdict(int)
    for read_name in read_names:
        read = read_dict[read_name]
        if read.alignments:
            read_lengths[read.get_length()] += 1
    estimated_genome_size = graph.get_estimated_sequence_len()

    # Now we need to finalise the bridges. This is the intensive step, as it involves creating a
    # consensus sequence, finding graph paths and doing alignments between the consensus and the
    # graph paths. We can therefore use threads to make this faster.
    long_read_bridges = [x for x in all_bridges if isinstance(x, LongReadBridge)]
    num_long_read_bridges = len(long_read_bridges)

    # We want to display this table one row at a time, so we have to fix all of the column widths
    # at the start.
    alignments, col_widths = get_long_read_bridge_table_parameters(graph, num_long_read_bridges,
                                                                   verbosity)
    if verbosity > 0:
        print_long_read_bridge_table_header(alignments, col_widths, verbosity)
    completed_count = 0
    if threads == 1:
        for bridge in long_read_bridges:
            output = bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                                     estimated_genome_size, expected_linear_seqs)
            completed_count += 1
            if verbosity > 0:
                print_long_read_bridge_table_row(alignments, col_widths, output, completed_count,
                                                 num_long_read_bridges, min_bridge_qual, verbosity)
    else:
        pool = ThreadPool(threads)
        arg_list = []

        # Sort the bridges based on how long they're predicted to take to finalise. This will make
        # the big ones runs first which helps to more efficiently use the CPU cores.
        # E.g. if the biggest bridge was at the end, we'd be left waiting for it to finish with
        # only one core (bad), but if it was at the start, other work could be done in parallel.
        long_read_bridges = sorted(long_read_bridges, reverse=True,
                                   key=lambda x: x.predicted_time_to_finalise())

        for bridge in long_read_bridges:
            arg_list.append((bridge, scoring_scheme, min_alignment_length, read_lengths,
                             estimated_genome_size, expected_linear_seqs))

        for output in pool.imap_unordered(finalise_bridge, arg_list):
            completed_count += 1
            if verbosity > 0:
                print_long_read_bridge_table_row(alignments, col_widths, output, completed_count,
                                                 num_long_read_bridges, min_bridge_qual, verbosity)
    return all_bridges


def get_long_read_bridge_table_parameters(graph, num_long_read_bridges, verbosity):
    max_seg_num_len = len(str(max(graph.segments.keys()))) + 1

    table_alignments = 'RL'
    table_col_widths = [2 * len(str(num_long_read_bridges)) + 1,  # Number
                        max_seg_num_len + 10]                     # Start to end

    if verbosity > 1:
        table_alignments += 'RR'
        table_col_widths += [5, 9]   # Read count, Consensus length

    if verbosity > 2:
        table_alignments += 'RR'
        table_col_widths += [9, 8]   # Consensus time, Target length

    if verbosity > 1:
        table_alignments += 'LRR'
        table_col_widths += [11, 8, 5]   # Search type, Search time, Path count

    table_alignments += 'L'
    table_col_widths += [40]  # Best path

    if verbosity > 2:
        table_alignments += 'RRRR'
        table_col_widths += [9, 9, 12, 11]  # Best path length, raw score, scaled score, len disc

    table_alignments += 'R'
    table_col_widths += [7]  # Quality

    return table_alignments, table_col_widths


def print_long_read_bridge_table_header(alignments, col_widths, verbosity):
    header_line_1 = ['', '']
    header_line_2 = ['', 'Start ' + get_right_arrow() + ' end']

    if verbosity > 1:
        header_line_1 += ['',      'Consensus']
        header_line_2 += ['Reads', 'len (bp)']
    if verbosity > 2:
        header_line_1 += ['Consensus', 'Target']
        header_line_2 += ['time (s)',  'len (bp)']
    if verbosity > 1:
        header_line_1 += ['',            'Search',   'Path']
        header_line_2 += ['Search type', 'time (s)', 'count']

    header_line_1 += ['']
    header_line_2 += ['Best path']

    if verbosity > 2:
        header_line_1 += ['Best path', 'Best path', 'Best path',    'Best path']
        header_line_2 += ['len (bp)',  'raw score', 'scaled score', 'length disc']

    header_line_1 += ['']
    header_line_2 += ['Quality']

    if any(x for x in header_line_1):
        print_table([header_line_1], col_separation=2, alignments=alignments,
                    header_format='normal', fixed_col_widths=col_widths, indent=0)
    print_table([header_line_2], col_separation=2, alignments=alignments,
                header_format='underline', fixed_col_widths=col_widths, indent=0)


def print_long_read_bridge_table_row(alignments, col_widths, output, completed_count,
                                     num_long_read_bridges, min_bridge_qual, verbosity):
    fraction = str(completed_count) + '/' + str(num_long_read_bridges)

    start, end, read_count, consensus_length, consensus_time, target_length, path_count, \
        search_type, search_time, best_path, best_path_len, best_path_raw_score, \
        best_path_scaled_score, best_path_length_discrepancy, quality = output

    start_to_end = (start + ' ' + get_right_arrow()).rjust(7) + ' ' + end
    quality_str = float_to_str(quality, 3)

    table_row = [fraction, start_to_end]

    if verbosity > 1:
        table_row += [read_count, consensus_length]

    if verbosity > 2:
        table_row += [consensus_time, target_length]

    if verbosity > 1:
        table_row += [search_type, search_time, path_count]

    table_row += [best_path]

    if verbosity > 2:
        table_row += [best_path_len, best_path_raw_score, best_path_scaled_score,
                      best_path_length_discrepancy]

    table_row += [quality_str]

    sub_colour = {}
    if quality < min_bridge_qual:
        sub_colour[quality_str] = 'red'
    print_table([table_row], col_separation=2, header_format='normal', indent=0,
                left_align_header=False, alignments=alignments, fixed_col_widths=col_widths,
                sub_colour=sub_colour)


def get_single_copy_alignments(read, single_copy_num_set, min_scaled_score):
    """
    Returns a list of single-copy segment alignments for the read.
    """
    sc_alignments = []
    for alignment in read.alignments:
        if alignment.ref.number in single_copy_num_set and \
                        alignment.scaled_score >= min_scaled_score:
            sc_alignments.append(alignment)
    return sc_alignments


def finalise_bridge(all_args):
    """
    Just a one-argument version of bridge.finalise, for pool.imap.
    """
    bridge, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size,\
        expected_linear_seqs = all_args
    return bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                           estimated_genome_size, expected_linear_seqs)


def reduce_expected_count(expected_count, a, b):
    """
    This function reduces the expected read count. It reduces by a factor which is a function of
    the read count, so low expected values aren't reduced much, but high expected values are
    reduced more. This is to help with high read depth cases where expected counts get quite high.

    https://www.desmos.com/calculator
    y=x\cdot \left(\left(\frac{a}{a+x}\right)\cdot \left(1-b\right)+b\right)
    """
    return expected_count * ((a / (a + expected_count)) * (1.0 - b) + b)
