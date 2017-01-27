"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Unicycler makes use of several C++ functions which are in cpp_functions.so. This module uses ctypes
to wrap them in similarly named Python functions.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
from ctypes import CDLL, cast, c_char_p, c_int, c_double, c_void_p, c_bool, POINTER
from .misc import quit_with_error


SO_FILE = 'cpp_functions.so'
SO_FILE_FULL = os.path.join(os.path.dirname(os.path.realpath(__file__)), SO_FILE)
if not os.path.isfile(SO_FILE_FULL):
    quit_with_error('could not find ' + SO_FILE + ' - please reinstall')
C_LIB = CDLL(SO_FILE_FULL)

# This is the big semi-global C++ Seqan alignment function at the heart of the aligner.
C_LIB.semiGlobalAlignment.argtypes = [c_char_p,  # Read name
                                      c_char_p,  # Read sequence
                                      c_int,     # Verbosity
                                      c_char_p,  # Minimap alignment info
                                      c_void_p,  # KmerPositions pointer
                                      c_int,     # Match score
                                      c_int,     # Mismatch score
                                      c_int,     # Gap open score
                                      c_int,     # Gap extension score
                                      c_double,  # Low score threshold
                                      c_bool,    # Return bad alignments
                                      c_int]     # Sensitivity level
C_LIB.semiGlobalAlignment.restype = c_void_p     # String describing alignments


def semi_global_alignment(read_name, read_sequence, verbosity, minimap_alignments_str,
                          kmer_positions_ptr, match_score, mismatch_score, gap_open_score,
                          gap_extend_score, low_score_threshold, keep_bad, sensitivity_level):
    """
    Python wrapper for semiGlobalAlignment C++ function.
    """
    ptr = C_LIB.semiGlobalAlignment(read_name.encode('utf-8'), read_sequence.encode('utf-8'),
                                    verbosity, minimap_alignments_str.encode('utf-8'),
                                    kmer_positions_ptr, match_score, mismatch_score,
                                    gap_open_score, gap_extend_score, low_score_threshold,
                                    keep_bad, sensitivity_level)
    return c_string_to_python_string(ptr)


# This is the global alignment function mainly used to compare read consensus sequences to assembly
# graph paths.
C_LIB.fullyGlobalAlignment.argtypes = [c_char_p,  # Sequence 1
                                       c_char_p,  # Sequence 2
                                       c_int,  # Match score
                                       c_int,  # Mismatch score
                                       c_int,  # Gap open score
                                       c_int,  # Gap extension score
                                       c_bool,  # Use banding
                                       c_int]  # Band size
C_LIB.fullyGlobalAlignment.restype = c_void_p  # String describing alignment


def fully_global_alignment(sequence_1, sequence_2, scoring_scheme, use_banding, band_size):
    """
    Python wrapper for fullyGlobalAlignment C++ function.
    """
    ptr = C_LIB.fullyGlobalAlignment(sequence_1.encode('utf-8'), sequence_2.encode('utf-8'),
                                     scoring_scheme.match, scoring_scheme.mismatch,
                                     scoring_scheme.gap_open, scoring_scheme.gap_extend,
                                     use_banding, band_size)
    return c_string_to_python_string(ptr)


# This is the mostly-global alignment function mainly used to compare potential path sequences to
# a read consensus. It is 'mostly-global' because there are free end gaps in the first sequence,
# so the path isn't penalised for not being complete.
C_LIB.pathAlignment.argtypes = [c_char_p,  # Sequence 1
                                c_char_p,  # Sequence 2
                                c_int,  # Match score
                                c_int,  # Mismatch score
                                c_int,  # Gap open score
                                c_int,  # Gap extension score
                                c_bool,  # Use banding
                                c_int]  # Band size
C_LIB.pathAlignment.restype = c_void_p  # String describing alignment


def path_alignment(partial_seq, full_seq, scoring_scheme, use_banding, band_size):
    """
    Python wrapper for pathAlignment C++ function.
    """
    ptr = C_LIB.pathAlignment(partial_seq.encode('utf-8'), full_seq.encode('utf-8'),
                              scoring_scheme.match, scoring_scheme.mismatch,
                              scoring_scheme.gap_open, scoring_scheme.gap_extend,
                              use_banding, band_size)
    return c_string_to_python_string(ptr)


# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None


def c_string_to_python_string(c_string):
    """
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    """
    python_string = cast(c_string, c_char_p).value.decode()
    C_LIB.freeCString(c_string)
    return python_string


# These functions make/delete a C++ object that will hold reference sequences for quick access.
C_LIB.newRefSeqs.argtypes = []
C_LIB.newRefSeqs.restype = c_void_p


def new_ref_seqs():
    """Python wrapper for newRefSeqs C++ function."""
    return C_LIB.newRefSeqs()


C_LIB.addRefSeq.argtypes = [c_void_p,  # SeqMap pointer
                            c_char_p,  # Name
                            c_char_p]  # Sequence
C_LIB.addRefSeq.restype = None


def add_ref_seq(ref_seqs_ptr, name, sequence):
    """Python wrapper for addRefSeq C++ function."""
    C_LIB.addRefSeq(ref_seqs_ptr, name.encode('utf-8'), sequence.encode('utf-8'))


C_LIB.deleteRefSeqs.argtypes = [c_void_p]
C_LIB.deleteRefSeqs.restype = None


def delete_ref_seqs(ref_seqs_ptr):
    """Python wrapper for deleteRefSeqs C++ function."""
    C_LIB.deleteRefSeqs(ref_seqs_ptr)


# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentScores.argtypes = [c_int,  # Random sequence length
                                                   c_int,  # Count
                                                   c_int,  # Match score
                                                   c_int,  # Mismatch score
                                                   c_int,  # Gap open score
                                                   c_int]  # Gap extension score
C_LIB.getRandomSequenceAlignmentScores.restype = c_void_p


def get_random_sequence_alignment_mean_and_std_dev(seq_length, count, scoring_scheme):
    """
    Python wrapper for getRandomSequenceAlignmentScores C++ function.
    """
    ptr = C_LIB.getRandomSequenceAlignmentScores(seq_length, count,
                                                 scoring_scheme.match, scoring_scheme.mismatch,
                                                 scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return_str = c_string_to_python_string(ptr)
    return_parts = return_str.split(',')
    return float(return_parts[0]), float(return_parts[1])


# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentErrorRates.argtypes = [c_int,  # Random sequence length
                                                       c_int,  # Count
                                                       c_int,  # Match score
                                                       c_int,  # Mismatch score
                                                       c_int,  # Gap open score
                                                       c_int]  # Gap extension score
C_LIB.getRandomSequenceAlignmentErrorRates.restype = c_void_p


def get_random_sequence_alignment_error_rates(seq_length, count, scoring_scheme):
    """
    Python wrapper for getRandomSequenceAlignmentErrorRate C++ function.
    """
    ptr = C_LIB.getRandomSequenceAlignmentErrorRates(seq_length, count,
                                                     scoring_scheme.match, scoring_scheme.mismatch,
                                                     scoring_scheme.gap_open,
                                                     scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)


# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.simulateDepths.argtypes = [POINTER(c_int),  # Alignment lengths
                                 c_int,  # Alignment count
                                 c_int,  # Reference length
                                 c_int,  # Iterations
                                 c_int]  # Threads
C_LIB.simulateDepths.restype = c_void_p


def simulate_depths(read_lengths, ref_length, iterations, threads):
    """
    Python wrapper for simulateDepths C++ function.
    """
    read_lengths_array = (c_int * len(read_lengths))(*read_lengths)
    ptr = C_LIB.simulateDepths(read_lengths_array, len(read_lengths), ref_length, iterations,
                               threads)
    return c_string_to_python_string(ptr)


# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.multipleSequenceAlignment.argtypes = [POINTER(c_char_p),  # Full-span sequences
                                            POINTER(c_char_p),  # Full-span qualities
                                            c_int,  # Full-span count
                                            POINTER(c_char_p),  # Start-only sequences
                                            POINTER(c_char_p),  # Start-only qualities
                                            c_int,  # Start-only count
                                            POINTER(c_char_p),  # End-only sequences
                                            POINTER(c_char_p),  # End-only qualities
                                            c_int,  # End-only count
                                            c_int,  # Bandwidth
                                            c_int,  # Match score
                                            c_int,  # Mismatch score
                                            c_int,  # Gap open score
                                            c_int]  # Gap extension score
C_LIB.multipleSequenceAlignment.restype = c_void_p


def multiple_sequence_alignment(full_length_sequences, full_length_qualities,
                                start_only_sequences, start_only_qualities,
                                end_only_sequences, end_only_qualities,
                                scoring_scheme, bandwidth=1000):
    """
    Python wrapper for multipleSequenceAlignment C++ function.
    """
    full_count = len(full_length_sequences)
    start_count = len(start_only_sequences)
    end_count = len(end_only_sequences)

    # At least one full length sequence is required.
    if not full_count:
        return "", [], [], []

    if len(full_length_qualities) < len(full_length_sequences):
        full_length_qualities += [""] * (len(full_length_sequences) - len(full_length_qualities))
    if len(start_only_qualities) < len(start_only_sequences):
        start_only_qualities += [""] * (len(start_only_sequences) - len(start_only_qualities))
    if len(end_only_qualities) < len(end_only_sequences):
        end_only_qualities += [""] * (len(end_only_sequences) - len(end_only_qualities))

    full_length_sequences = [x.encode('utf-8') for x in full_length_sequences]
    full_length_qualities = [x.encode('utf-8') for x in full_length_qualities]
    start_only_sequences = [x.encode('utf-8') for x in start_only_sequences]
    start_only_qualities = [x.encode('utf-8') for x in start_only_qualities]
    end_only_sequences = [x.encode('utf-8') for x in end_only_sequences]
    end_only_qualities = [x.encode('utf-8') for x in end_only_qualities]

    sequences_1 = (c_char_p * len(full_length_sequences))(*full_length_sequences)
    qualities_1 = (c_char_p * len(full_length_qualities))(*full_length_qualities)
    sequences_2 = (c_char_p * len(start_only_sequences))(*start_only_sequences)
    qualities_2 = (c_char_p * len(start_only_qualities))(*start_only_qualities)
    sequences_3 = (c_char_p * len(end_only_sequences))(*end_only_sequences)
    qualities_3 = (c_char_p * len(end_only_qualities))(*end_only_qualities)

    ptr = C_LIB.multipleSequenceAlignment(sequences_1, qualities_1, full_count,
                                          sequences_2, qualities_2, start_count,
                                          sequences_3, qualities_3, end_count,
                                          bandwidth,
                                          scoring_scheme.match, scoring_scheme.mismatch,
                                          scoring_scheme.gap_open, scoring_scheme.gap_extend)
    result = c_string_to_python_string(ptr)
    result_parts = result.split(';')

    consensus = result_parts[0]

    all_scores = result_parts[1].split(',')
    full_length_scores = all_scores[:full_count]
    start_only_scores = all_scores[full_count:full_count + start_count]
    if end_count:
        end_only_scores = all_scores[-end_count:]
    else:
        end_only_scores = []

    # between_seq_scores = [x.split(',') for x in result_parts[2:]]

    return consensus, full_length_scores, start_only_scores, end_only_scores


# This function conducts a minimap alignment between reads and reference.
C_LIB.minimapAlignReads.argtypes = [c_char_p,  # Reference fasta filename
                                    c_char_p,  # Reads fastq filename
                                    c_int,     # Threads
                                    c_int,     # Sensitivity level
                                    c_int]     # Verbosity
C_LIB.minimapAlignReads.restype = c_void_p     # String describing alignments


def minimap_align_reads(reference_fasta, reads_fastq, threads, verbosity, sensitivity_level):
    """
    Python wrapper for destroyMinimapIndex C++ function.
    """
    ptr = C_LIB.minimapAlignReads(reference_fasta.encode('utf-8'), reads_fastq.encode('utf-8'),
                                  threads, verbosity, sensitivity_level)
    return c_string_to_python_string(ptr)
