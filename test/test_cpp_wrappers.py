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
import unicycler.cpp_wrappers
import unicycler.read_ref
import unicycler.alignment
import unicycler.misc


class TestFullyGlobalAlignment(unittest.TestCase):

    def setUp(self):
        test_fasta = os.path.join(os.path.dirname(__file__), 'test_1.fasta')
        fasta = unicycler.misc.load_fasta(test_fasta)
        self.seqs = [x[1] for x in fasta]
        self.scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')

    @staticmethod
    def get_raw_and_scaled_scores(result):
        seqan_parts = result.split(',', 9)
        raw_score = int(seqan_parts[6])
        scaled_score = float(seqan_parts[7])
        return raw_score, scaled_score

    def test_perfect_alignment(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[1],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 60)
        self.assertEqual(scaled_score, 100.0)

    def test_one_mismatch(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[2],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 51)
        self.assertTrue(scaled_score < 100.0)

    def test_1bp_insertion(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[3],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 55)
        self.assertTrue(scaled_score < 100.0)

    def test_1bp_deletion(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[4],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 52)
        self.assertTrue(scaled_score < 100.0)

    def test_2bp_insertion(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[5],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 53)
        self.assertTrue(scaled_score < 100.0)

    def test_2bp_deletion(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[6],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 47)
        self.assertTrue(scaled_score < 100.0)

    def test_2bp_insertion_and_2bp_deletion(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[0], self.seqs[7],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 40)
        self.assertTrue(scaled_score < 100.0)

    def test_long_perfect_alignment(self):
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[8], self.seqs[9],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 3000)
        self.assertEqual(scaled_score, 100.0)

    def test_shift_big_band(self):
        """
        This test has a 20 bp insertion near the start and a 20 bp deletion near the end, shifting
        the alignment by 20 bp in the middle.
        """
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[8], self.seqs[10],
                                                                self.scoring_scheme, True, 1000)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertEqual(raw_score, 2854)

    def test_shift_small_band(self):
        """
        This test uses the same sequences as the previous test, but uses a smaller bandwidth than
        the sequence shift induced by the indels, therefore producing a worse alignment.
        """
        result = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[8], self.seqs[10],
                                                                self.scoring_scheme, True, 10)
        raw_score, scaled_score = self.get_raw_and_scaled_scores(result)
        self.assertTrue(raw_score < 2854)

    def test_random_seqs(self):
        """
        This test aligns two random sequences of different length, making sure that the scores are
        the same regardless of which comes first.
        """
        result_1 = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[11], self.seqs[12],
                                                                  self.scoring_scheme, True, 1000)
        raw_score_1, scaled_score_1 = self.get_raw_and_scaled_scores(result_1)
        result_2 = unicycler.cpp_wrappers.fully_global_alignment(self.seqs[12], self.seqs[11],
                                                                  self.scoring_scheme, True, 1000)
        raw_score_2, scaled_score_2 = self.get_raw_and_scaled_scores(result_2)
        self.assertEqual(raw_score_1, raw_score_2)
        self.assertEqual(scaled_score_1, scaled_score_2)


class TestPathAlignment(unittest.TestCase):
    pass


class TestMultipleSequenceAlignment(unittest.TestCase):

    def setUp(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_2.fastq')
        self.read_dict, self.read_names, _ = unicycler.read_ref.load_long_reads(test_fastq, 0)
        self.scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        self.seqs = [self.read_dict[name].sequence for name in self.read_names]
        self.quals = [self.read_dict[name].qualities for name in self.read_names]
        self.original_seq = self.seqs[0]

    def test_consensus_with_subs(self):
        seqs = self.seqs[1:4]
        quals = self.quals[1:4]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_deletions(self):
        seqs = self.seqs[4:7]
        quals = self.quals[4:7]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_insertions(self):
        seqs = self.seqs[7:10]
        quals = self.quals[7:10]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_deletions_and_insertions(self):
        seqs = self.seqs[4:10]
        quals = self.quals[4:10]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_all(self):
        seqs = self.seqs[1:10]
        quals = self.quals[1:10]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_two_way_consensus(self):
        seqs = self.seqs[10:12]
        quals = self.quals[10:12]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_different_qualities(self):
        seqs = self.seqs[12:16]
        quals = self.quals[12:16]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)
        self.assertEqual(scores[0], 1.0)
        self.assertTrue(scores[0] > scores[1])
        self.assertTrue(scores[1] > scores[2])
        self.assertTrue(scores[2] > scores[3])

    def test_start_end_insertions(self):
        seqs = [self.seqs[0]] + self.seqs[16:18]
        quals = [self.quals[0]] + self.quals[16:18]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_start_end_deletions(self):
        seqs = [self.seqs[0]] + self.seqs[18:20]
        quals = [self.quals[0]] + self.quals[18:20]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_start_end_insertions_and_deletions(self):
        seqs = self.seqs[16:20]
        quals = self.quals[16:20]
        consensus, scores = unicycler.cpp_wrappers.consensus_alignment(seqs, quals,
                                                                        self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)
