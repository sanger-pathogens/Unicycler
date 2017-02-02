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
        self.verbosity = 0
        temp_name = 'TEMP_' + str(os.getpid())
        self.temp_fasta = temp_name + '.fasta'
        self.temp_fastq = temp_name + '.fastq'
        all_ref_fasta = os.path.join(os.path.dirname(__file__),
                                     'test_semi_global_alignment_tough.fasta')
        self.all_refs = unicycler.read_ref.load_references(all_ref_fasta, self.verbosity)
        all_read_fastq = os.path.join(os.path.dirname(__file__),
                                           'test_semi_global_alignment_tough.fastq')
        self.all_reads, _, _ = unicycler.read_ref.load_long_reads(all_read_fastq, self.verbosity)

        # Alignments are tests with approximate boundaries to allow for a big of wiggle room if the
        # implementation changes.
        self.pos_margin_of_error = 20

    def do_alignment(self, read_ref_name):
        # Save just the read/ref of interest (which will have the same name) into temporary
        # separate files, so th alignment can be done for just them.
        ref = [x for x in self.all_refs if x.name == read_ref_name][0]
        read = [x for x in self.all_reads.values() if x.name == read_ref_name][0]
        with open(self.temp_fasta, 'wt') as ref_fasta:
            ref_fasta.write('>' + ref.name + '\n')
            ref_fasta.write(ref.sequence + '\n')
        with open(self.temp_fastq, 'wt') as read_fastq:
            read_fastq.write('@' + read.name + '\n')
            read_fastq.write(read.sequence + '\n')
            read_fastq.write('+\n')
            read_fastq.write(read.qualities + '\n')

        refs = unicycler.read_ref.load_references(self.temp_fasta, self.verbosity)
        read_dict, read_names, _ = unicycler.read_ref.load_long_reads(self.temp_fastq,
                                                                      self.verbosity)
        scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        sensitivity_level = 0
        contamination_fasta = None
        threads = 1
        min_align_length = 10
        allowed_overlap = 0
        self.aligned_reads = unicycler.unicycler_align.\
                semi_global_align_long_reads(refs, self.temp_fasta, read_dict, read_names,
                                             self.temp_fastq, threads, scoring_scheme, [None],
                                             False, min_align_length, None, None, allowed_overlap,
                                             sensitivity_level, contamination_fasta, self.verbosity)

    def tearDown(self):
        if os.path.isfile(self.temp_fasta):
            os.remove(self.temp_fasta)
        if os.path.isfile(self.temp_fastq):
            os.remove(self.temp_fastq)

    def test_tough_alignment_0(self):
        """
        The beginning of the reference in this case is repetitive, which was able to throw off
        Seqan's global chaining algorithm, resulting in an awkward alignment. I think I fixed this
        by limiting Seqan to the seeds which are near the diagonals of line tracing points.
        """
        self.do_alignment('0')
        read = self.aligned_reads['0']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '0')
        self.assertTrue(alignment.raw_score >= 126525)
        self.assertTrue(alignment.scaled_score > 91.19)
        read_start, read_end = alignment.read_start_end_positive_strand()
        self.assertTrue(abs(read_start - 18662) < self.pos_margin_of_error)
        self.assertEqual(read_end, 72402)  # end of read
        self.assertEqual(alignment.ref_start_pos, 0)     # start of ref
        self.assertTrue(abs(alignment.ref_end_pos - 55814) < self.pos_margin_of_error)

    def test_tough_alignment_1(self):
        """
        This read goes through a repetitive area, which means that the densest area of common
        k-mers is not on the correct alignment line. This means more than one line tracing is
        required to get it right.
        """
        self.do_alignment('1')
        read = self.aligned_reads['1']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '1')
        self.assertTrue(alignment.raw_score >= 20740)
        self.assertTrue(alignment.scaled_score > 91.02)
        read_start, read_end = alignment.read_start_end_positive_strand()
        self.assertTrue(abs(read_start - 10785) < self.pos_margin_of_error)
        self.assertTrue(abs(read_end - 19629) < self.pos_margin_of_error)
        self.assertEqual(alignment.ref_start_pos, 0)   # start of ref
        self.assertEqual(alignment.ref_end_pos, 9241)  # end of ref

    def test_tough_alignment_2(self):
        """
        This read goes through a repetitive area, which means that the densest area of common
        k-mers is not on the correct alignment line. This means more than one line tracing is
        required to get it right.
        """
        self.do_alignment('2')
        read = self.aligned_reads['2']
        self.assertEqual(len(read.alignments), 1)
        alignment = read.alignments[0]
        self.assertEqual(alignment.read.name, '2')
        self.assertTrue(alignment.raw_score >= 34449)
        self.assertTrue(alignment.scaled_score > 90.35)
        read_start, read_end = alignment.read_start_end_positive_strand()
        self.assertTrue(abs(read_start - 22493) < self.pos_margin_of_error)
        self.assertEqual(read_end, 37581)  # end of read
        self.assertEqual(alignment.ref_start_pos, 0)     # start of ref
        self.assertTrue(abs(alignment.ref_end_pos - 15673) < self.pos_margin_of_error)
