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
import shutil
import unicycler.misc
import unicycler.blast_func
import unicycler.assembly_graph_segment


class TestBlastFunc(unittest.TestCase):
    """
    Tests rotation of completed replicons to a predefined start.
    """

    def setUp(self):
        test_fasta = os.path.join(os.path.dirname(__file__), 'test_blast_func_sequences.fasta')
        self.fasta = unicycler.misc.load_fasta(test_fasta)
        self.start_genes = os.path.join(os.path.dirname(__file__),
                                        'test_blast_func_start_genes.fasta')
        self.threads = 1
        self.verbosity = 0
        self.start_gene_id = 90.0
        self.start_gene_cov = 95.0
        self.blast_dir = 'TEMP_' + str(os.getpid())
        if not os.path.exists(self.blast_dir):
            os.makedirs(self.blast_dir)

    def tearDown(self):
        if os.path.exists(self.blast_dir):
            shutil.rmtree(self.blast_dir)

    def test_random_seq_no_start_gene(self):
        seq = [x for x in self.fasta if x[0] == 'random_seq_no_start_gene'][0][1]
        with self.assertRaises(unicycler.blast_func.CannotFindStart):
            unicycler.blast_func.find_start_gene(seq, self.start_genes, self.start_gene_id,
                                                 self.start_gene_cov, self.blast_dir,
                                                 'makeblastdb', 'tblastn', self.threads)

    def test_random_seq_with_exact_gene_forward_strand(self):
        seq = [x for x in self.fasta if x[0] == 'random_seq_with_exact_gene_forward_strand'][0][1]
        hit = unicycler.blast_func.find_start_gene(seq, self.start_genes, self.start_gene_id,
                                                   self.start_gene_cov, self.blast_dir,
                                                   'makeblastdb', 'tblastn', self.threads)
        self.assertEqual(hit.qseqid, 'UniRef90_P66818')
        self.assertEqual(hit.start_pos, 36661)
        self.assertFalse(hit.flip)
        self.assertEqual(hit.pident, 100.0)
        self.assertEqual(hit.query_cov, 100.0)

        length_before_rotate = len(seq)
        seg = unicycler.assembly_graph_segment.Segment(1, 1.0, seq, True)
        seg.rotate_sequence(hit.start_pos, hit.flip)

        self.assertEqual(len(seg.forward_sequence), length_before_rotate)
        self.assertTrue(seg.forward_sequence.startswith('ATGCAGGAACGCATTAAAGCGTGCTTTACCGAAAG'))

    def test_random_seq_with_exact_gene_reverse_strand(self):
        seq = [x for x in self.fasta if x[0] == 'random_seq_with_exact_gene_reverse_strand'][0][1]
        hit = unicycler.blast_func.find_start_gene(seq, self.start_genes, self.start_gene_id,
                                                   self.start_gene_cov, self.blast_dir,
                                                   'makeblastdb', 'tblastn', self.threads)
        self.assertEqual(hit.qseqid, 'UniRef90_P66818')
        self.assertEqual(hit.start_pos, 82415)
        self.assertTrue(hit.flip)
        self.assertEqual(hit.pident, 100.0)
        self.assertEqual(hit.query_cov, 100.0)

        length_before_rotate = len(seq)
        seg = unicycler.assembly_graph_segment.Segment(1, 1.0, seq, True)
        seg.rotate_sequence(hit.start_pos, hit.flip)

        self.assertEqual(len(seg.forward_sequence), length_before_rotate)
        self.assertTrue(seg.forward_sequence.startswith('ATGCAGGAACGCATTAAAGCGTGCTTTACCGAAAG'))
