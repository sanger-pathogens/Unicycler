"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import unicycler.read_ref
import unicycler.alignment
import unicycler.unicycler_align


class TestPerfectMatchAlignments(unittest.TestCase):

    def setUp(self):
        ref_fasta = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fasta')
        read_fastq = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fastq')
        verbosity = 0
        refs = unicycler.read_ref.load_references(ref_fasta, verbosity)
        read_dict, read_names, _ = unicycler.read_ref.load_long_reads(read_fastq, verbosity)
        scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        sensitivity_level = 0
        contamination_fasta = None
        threads = 1
        min_align_length = 10
        allowed_overlap = 0
        self.aligned_reads = unicycler.unicycler_align.\
                semi_global_align_long_reads(refs, ref_fasta, read_dict, read_names, read_fastq,
                                             threads, scoring_scheme, [None], False,
                                             min_align_length, None, None, allowed_overlap,
                                             sensitivity_level, contamination_fasta, verbosity)

    def test_read_contained_1(self):
        read = self.aligned_reads['0']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '0')
        self.assertEqual(alignment.raw_score, 300)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 100)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 100)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 60)
        self.assertEqual(alignment.ref_end_pos, 160)
        self.assertEqual(len(alignment.cigar_parts), 1)
        self.assertEqual(alignment.cigar_parts[0], '100M')

    def test_read_contained_2(self):
        read = self.aligned_reads['1']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '1')
        self.assertEqual(alignment.raw_score, 600)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 200)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 200)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 100)
        self.assertEqual(alignment.ref_end_pos, 300)
        self.assertEqual(len(alignment.cigar_parts), 1)
        self.assertEqual(alignment.cigar_parts[0], '200M')

    def test_read_contained_3(self):
        read = self.aligned_reads['2']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '2')
        self.assertEqual(alignment.raw_score, 450)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 150)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 150)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 150)
        self.assertEqual(len(alignment.cigar_parts), 1)
        self.assertEqual(alignment.cigar_parts[0], '150M')

    def test_ref_contained_1(self):
        read = self.aligned_reads['3']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '3')
        self.assertEqual(alignment.raw_score, 300)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 100)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 62)
        self.assertEqual(alignment.read_end_pos, 162)
        self.assertEqual(alignment.read_end_gap, 138)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 100)
        self.assertEqual(len(alignment.cigar_parts), 3)
        self.assertEqual(alignment.cigar_parts[0], '62S')
        self.assertEqual(alignment.cigar_parts[1], '100M')
        self.assertEqual(alignment.cigar_parts[2], '138S')

    def test_ref_contained_2(self):
        read = self.aligned_reads['4']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '4')
        self.assertEqual(alignment.raw_score, 360)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 120)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 120)
        self.assertEqual(alignment.read_end_gap, 180)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 120)
        self.assertEqual(len(alignment.cigar_parts), 2)
        self.assertEqual(alignment.cigar_parts[0], '120M')
        self.assertEqual(alignment.cigar_parts[1], '180S')

    def test_ref_contained_3(self):
        read = self.aligned_reads['5']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '5')
        self.assertEqual(alignment.raw_score, 540)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 180)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 120)
        self.assertEqual(alignment.read_end_pos, 300)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 180)
        self.assertEqual(len(alignment.cigar_parts), 2)
        self.assertEqual(alignment.cigar_parts[0], '120S')
        self.assertEqual(alignment.cigar_parts[1], '180M')

    def test_read_start_overlap(self):
        read = self.aligned_reads['6']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '6')
        self.assertEqual(alignment.raw_score, 330)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 110)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 190)
        self.assertEqual(alignment.read_end_pos, 300)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 110)
        self.assertEqual(len(alignment.cigar_parts), 2)
        self.assertEqual(alignment.cigar_parts[0], '190S')
        self.assertEqual(alignment.cigar_parts[1], '110M')

    def test_read_end_overlap(self):
        read = self.aligned_reads['7']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '7')
        self.assertEqual(alignment.raw_score, 390)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 130)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 130)
        self.assertEqual(alignment.read_end_gap, 170)
        self.assertEqual(alignment.ref_start_pos, 170)
        self.assertEqual(alignment.ref_end_pos, 300)
        self.assertEqual(len(alignment.cigar_parts), 2)
        self.assertEqual(alignment.cigar_parts[0], '130M')
        self.assertEqual(alignment.cigar_parts[1], '170S')

    def test_end_to_end(self):
        read = self.aligned_reads['8']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '8')
        self.assertEqual(alignment.raw_score, 900)
        self.assertEqual(alignment.scaled_score, 100.0)
        self.assertEqual(alignment.percent_identity, 100.0)
        self.assertEqual(alignment.match_count, 300)
        self.assertEqual(alignment.mismatch_count, 0)
        self.assertEqual(alignment.insertion_count, 0)
        self.assertEqual(alignment.deletion_count, 0)
        self.assertEqual(alignment.read_start_pos, 0)
        self.assertEqual(alignment.read_end_pos, 300)
        self.assertEqual(alignment.read_end_gap, 0)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertEqual(alignment.ref_end_pos, 300)
        self.assertEqual(len(alignment.cigar_parts), 1)
        self.assertEqual(alignment.cigar_parts[0], '300M')


class TestToughAlignments(unittest.TestCase):
    """
    These test cases are made from real alignments which proved to be difficult.
    """

    def setUp(self):
        ref_fasta = os.path.join(os.path.dirname(__file__),
                                 'test_semi_global_alignment_tough.fasta')
        read_fastq = os.path.join(os.path.dirname(__file__),
                                  'test_semi_global_alignment_tough.fastq')
        verbosity = 0
        refs = unicycler.read_ref.load_references(ref_fasta, verbosity)
        read_dict, read_names, _ = unicycler.read_ref.load_long_reads(read_fastq, verbosity)
        scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        sensitivity_level = 0
        contamination_fasta = None
        threads = 1
        min_align_length = 10
        allowed_overlap = 0
        self.aligned_reads = unicycler.unicycler_align.\
                semi_global_align_long_reads(refs, ref_fasta, read_dict, read_names, read_fastq,
                                             threads, scoring_scheme, [None], False,
                                             min_align_length, None, None, allowed_overlap,
                                             sensitivity_level, contamination_fasta, verbosity)

    def test_tough_alignment_1(self):
        """
        The beginning of the reference in this case is repetitive, which was able to throw off
        Seqan's global chaining algorithm, resulting in an awkward alignment. I think I fixed this
        by limiting Seqan to the seeds which are near the diagonals of line tracing points.
        """
        read = self.aligned_reads['0']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '0')
        self.assertTrue(alignment.raw_score >= 126525)
        self.assertTrue(alignment.scaled_score > 91.19)
        self.assertTrue(abs(alignment.read_start_pos - 18662) < 20)
        self.assertEqual(alignment.read_end_pos, 72402)
        self.assertEqual(alignment.ref_start_pos, 0)
        self.assertTrue(abs(alignment.ref_end_pos - 55814) < 20)
