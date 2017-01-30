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
import unicycler.misc


class TestMiscFunctions(unittest.TestCase):

    def test_float_to_str(self):
        self.assertEqual('14.3', unicycler.misc.float_to_str(14.312, 1))
        self.assertEqual('1,000,000.1', unicycler.misc.float_to_str(1000000.101, 1))
        self.assertEqual('14.31', unicycler.misc.float_to_str(14.312, 2))
        self.assertEqual('14.32', unicycler.misc.float_to_str(14.319, 2))
        self.assertEqual(' 14.32', unicycler.misc.float_to_str(14.319, 2, 100))
        self.assertEqual(' 14.32', unicycler.misc.float_to_str(14.319, 2, 100.99999999))
        self.assertEqual('   14.32', unicycler.misc.float_to_str(14.319, 2, 999.9999999))

    def test_int_to_str(self):
        self.assertEqual('14', unicycler.misc.int_to_str(14))
        self.assertEqual('140', unicycler.misc.int_to_str(140))
        self.assertEqual('1,400', unicycler.misc.int_to_str(1400))
        self.assertEqual('14,000', unicycler.misc.int_to_str(14000))
        self.assertEqual('       14,000', unicycler.misc.int_to_str(14000, 1000000000))

    def test_get_nice_header(self):
        self.assertEqual(unicycler.misc.get_nice_header('NODE_1_length_49055_cov_41.580185'),
                         'NODE_1')
        self.assertEqual(unicycler.misc.get_nice_header('NODE_69_length_12837_cov_39.210640'),
                         'NODE_69')
        self.assertEqual(unicycler.misc.get_nice_header('name stuff'),
                         'name')
        self.assertEqual(unicycler.misc.get_nice_header('name stuff stuff stuff'),
                         'name')

    def test_is_header_spades_format(self):
        self.assertTrue(unicycler.misc.is_header_spades_format('NODE_1_length_49055_'
                                                               'cov_41.580185'))
        self.assertTrue(unicycler.misc.is_header_spades_format('NODE_69_length_12837_'
                                                               'cov_39.210640'))
        self.assertFalse(unicycler.misc.is_header_spades_format('name stuff'))
        self.assertFalse(unicycler.misc.is_header_spades_format('name stuff stuff stuff'))

    def test_reverse_complement(self):
        self.assertEqual('GCAGGCCGCTTAATGAATAGATCATGGCTGCGCCGCCTACCGGTCCGAGACCTTCGCTGA',
                         unicycler.misc.reverse_complement('TCAGCGAAGGTCTCGGACCGGTAGGCGGCGCAGCCATG'
                                                           'ATCTATTCATTAAGCGGCCTGC'))
        self.assertEqual('', unicycler.misc.reverse_complement(''))
        self.assertEqual('TATTTNGTTANAT', unicycler.misc.reverse_complement('ATNTAACNAAATA'))
        self.assertEqual('ATNTAACNAAATA', unicycler.misc.reverse_complement('TATTTNGTTANAT'))
        self.assertEqual('TGACBWDARAYACHASKGVTMACNG',
                         unicycler.misc.reverse_complement('CNGTKABCMSTDGTRTYTHWVGTCA'))
        self.assertEqual('tgACBWDARAYACHASKGVTMACnG',
                         unicycler.misc.reverse_complement('CnGTKABCMSTDGTRTYTHWVGTca'))

    def test_get_random_base(self):
        a_count, c_count, g_count, t_count, other_count = 0, 0, 0, 0, 0
        for i in range(10000):
            base = unicycler.misc.get_random_base()
            if base == 'A':
                a_count += 1
            elif base == 'C':
                c_count += 1
            elif base == 'G':
                g_count += 1
            elif base == 'T':
                t_count += 1
            else:
                other_count += 1
        self.assertTrue(a_count > 0)
        self.assertTrue(c_count > 0)
        self.assertTrue(g_count > 0)
        self.assertTrue(t_count > 0)
        self.assertTrue(other_count == 0)

    def test_get_random_sequence(self):
        self.assertEqual(10, len(unicycler.misc.get_random_sequence(10)))
        self.assertEqual(100, len(unicycler.misc.get_random_sequence(100)))
        self.assertEqual(1000, len(unicycler.misc.get_random_sequence(1000)))

    def test_get_percentile(self):
        self.assertEqual(20, unicycler.misc.get_percentile([50, 20, 40, 35, 15], 30))
        self.assertEqual(20, unicycler.misc.get_percentile([20, 50, 40, 35, 15], 40))
        self.assertEqual(35, unicycler.misc.get_percentile([50, 20, 40, 35, 15], 50))
        self.assertEqual(50, unicycler.misc.get_percentile([50, 20, 15, 35, 40], 100))
        self.assertEqual(7, unicycler.misc.get_percentile([3, 16, 7, 8, 8, 13, 10, 15, 6, 20],
                                                          25))
        self.assertEqual(8, unicycler.misc.get_percentile([16, 7, 8, 8, 13, 10, 15, 6, 20, 3],
                                                          50))
        self.assertEqual(15, unicycler.misc.get_percentile([3, 16, 7, 15, 8, 13, 10, 8, 6, 20],
                                                           75))
        self.assertEqual(20, unicycler.misc.get_percentile([20, 16, 7, 8, 8, 13, 10, 15, 6, 3],
                                                           100))
        self.assertEqual(7, unicycler.misc.get_percentile([3, 6, 7, 8, 8, 9, 10, 13, 15, 16, 20],
                                                          25))
        self.assertEqual(9, unicycler.misc.get_percentile([7, 9, 10, 3, 8, 15, 16, 13, 8, 20, 6],
                                                          50))
        self.assertEqual(15, unicycler.misc.get_percentile([3, 15, 13, 7, 8, 20, 10, 8, 16, 6, 9],
                                                           75))
        self.assertEqual(20, unicycler.misc.get_percentile([6, 13, 10, 8, 15, 3, 9, 8, 16, 7, 20],
                                                           100))

    def test_weighted_average(self):
        self.assertAlmostEqual(1.0, unicycler.misc.weighted_average(1.0, 2.0, 1.0, 0.0))
        self.assertAlmostEqual(2.0, unicycler.misc.weighted_average(1.0, 2.0, 0.0, 1.0))
        self.assertAlmostEqual(1.5, unicycler.misc.weighted_average(1.0, 2.0, 0.5, 0.5))
        self.assertAlmostEqual(1.5, unicycler.misc.weighted_average(1.0, 2.0, 1.5, 1.5))
        self.assertAlmostEqual(1.5, unicycler.misc.weighted_average(1.0, 2.0, 1.5, 1.5))
        self.assertAlmostEqual(1.3333333333333, unicycler.misc.weighted_average(1.0, 2.0, 4, 2))
        self.assertAlmostEqual(1.6666666666666, unicycler.misc.weighted_average(1.0, 2.0, 3, 6))
        self.assertAlmostEqual(1.6666666666666, unicycler.misc.weighted_average(1.0, 2.0, 3, 6))

    def test_weighted_average_list(self):
        self.assertAlmostEqual(1.0, unicycler.misc.weighted_average_list([1.0, 2.0, 3.0],
                                                                         [2.0, 0.0, 0.0]))
        self.assertAlmostEqual(2.0, unicycler.misc.weighted_average_list([1.0, 2.0, 3.0],
                                                                         [2.0, 0.0, 2.0]))
        self.assertAlmostEqual(2.5, unicycler.misc.weighted_average_list([1.0, 2.0, 3.0],
                                                                         [0.0, 1.0, 1.0]))
        self.assertAlmostEqual(2.6666666666666,
                               unicycler.misc.weighted_average_list([1.0, 2.0, 3.0],
                                                                    [0.0, 1.0, 2.0]))

    def test_round_to_nearest_odd(self):
        self.assertEqual(1, unicycler.misc.round_to_nearest_odd(0.9))
        self.assertEqual(1, unicycler.misc.round_to_nearest_odd(0.1))
        self.assertEqual(1, unicycler.misc.round_to_nearest_odd(1.9))
        self.assertEqual(3, unicycler.misc.round_to_nearest_odd(2.1))
        self.assertEqual(3, unicycler.misc.round_to_nearest_odd(2.9))
        self.assertEqual(3, unicycler.misc.round_to_nearest_odd(3.9))
        self.assertEqual(5, unicycler.misc.round_to_nearest_odd(4.1))
        self.assertEqual(5, unicycler.misc.round_to_nearest_odd(5.0))

    def test_get_compression_type(self):
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sample_data')
        ref = os.path.join(sample_dir, 'reference.fasta')
        reads = os.path.join(sample_dir, 'short_reads_1.fastq.gz')
        self.assertEqual('plain', unicycler.misc.get_compression_type(ref))
        self.assertEqual('gz', unicycler.misc.get_compression_type(reads))

    def test_get_sequence_file_type(self):
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sample_data')
        ref = os.path.join(sample_dir, 'reference.fasta')
        reads = os.path.join(sample_dir, 'short_reads_1.fastq.gz')
        self.assertEqual('FASTA', unicycler.misc.get_sequence_file_type(ref))
        self.assertEqual('FASTQ', unicycler.misc.get_sequence_file_type(reads))

    def test_get_num_agreement(self):
        self.assertAlmostEqual(1.0, unicycler.misc.get_num_agreement(2.0, 2.0))
        self.assertAlmostEqual(1.0, unicycler.misc.get_num_agreement(200.0, 200.0))
        self.assertAlmostEqual(0.0, unicycler.misc.get_num_agreement(0.0, 2.0))
        self.assertAlmostEqual(0.0, unicycler.misc.get_num_agreement(0.0, 200.0))
        self.assertAlmostEqual(0.5, unicycler.misc.get_num_agreement(1.0, 2.0))
        self.assertAlmostEqual(0.5, unicycler.misc.get_num_agreement(100.0, 200.0))

    def test_flip_number_order(self):
        self.assertEqual(unicycler.misc.flip_number_order(5, 6), ((5, 6), False))
        self.assertEqual(unicycler.misc.flip_number_order(6, 5), ((6, 5), False))
        self.assertEqual(unicycler.misc.flip_number_order(-5, -6), ((6, 5), True))
        self.assertEqual(unicycler.misc.flip_number_order(-6, -5), ((5, 6), True))

    def test_load_fasta(self):
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sample_data')
        ref = os.path.join(sample_dir, 'reference.fasta')
        fasta = unicycler.misc.load_fasta(ref)
        self.assertEqual(len(fasta), 3)
        self.assertEqual(fasta[0][0], 'NC_016833.1')
        self.assertEqual(fasta[2][0], 'NC_016834.1')
        self.assertTrue(fasta[0][1].startswith('ATGCTGATGAAAATACCTAAATAATCAGCCAGCACTCTATCTTTCCAAAT'
                                               'CCACAGCATAGCAAAGAGAGCAAAAGAGCCTGTAAATTCAGAAATT'))
        self.assertTrue(fasta[2][1].endswith('AGTTGATTTAAATCGCTACACCATTATGATTCATGTAGCGATTTAAATTACT'
                                             'ACATAATGGTGATTAGC'))

    def test_load_fasta_with_full_header(self):
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sample_data')
        ref = os.path.join(sample_dir, 'reference.fasta')
        fasta = unicycler.misc.load_fasta_with_full_header(ref)
        self.assertEqual(len(fasta), 3)
        self.assertEqual(fasta[0][0], 'NC_016833.1')
        self.assertEqual(fasta[2][0], 'NC_016834.1')
        self.assertEqual(fasta[0][1], 'NC_016833.1 Shigella sonnei 53G plasmid A, complete genome')
        self.assertEqual(fasta[2][1], 'NC_016834.1 Shigella sonnei 53G plasmid E, complete genome')
        self.assertTrue(fasta[0][2].startswith('ATGCTGATGAAAATACCTAAATAATCAGCCAGCACTCTATCTTTCCAAAT'
                                               'CCACAGCATAGCAAAGAGAGCAAAAGAGCCTGTAAATTCAGAAATT'))
        self.assertTrue(fasta[2][2].endswith('AGTTGATTTAAATCGCTACACCATTATGATTCATGTAGCGATTTAAATTACT'
                                             'ACATAATGGTGATTAGC'))

    def test_score_function(self):
        self.assertAlmostEqual(unicycler.misc.score_function(0.0, 1.0), 0.0)
        self.assertAlmostEqual(unicycler.misc.score_function(0.0, 2.0), 0.0)
        self.assertAlmostEqual(unicycler.misc.score_function(0.0, 100.0), 0.0)
        self.assertAlmostEqual(unicycler.misc.score_function(2.0, 2.0), 0.5)
        self.assertAlmostEqual(unicycler.misc.score_function(12.0, 12.0), 0.5)
        self.assertAlmostEqual(unicycler.misc.score_function(24.0, 12.0), 0.666666666666666666666)
        self.assertAlmostEqual(unicycler.misc.score_function(6.0, 12.0), 0.3333333333333333333333)
        self.assertAlmostEqual(unicycler.misc.score_function(100000000000000.0, 1.0), 1.0)
        self.assertAlmostEqual(unicycler.misc.score_function(100000000000000.0, 10.0), 1.0)

    def test_strip_read_extensions(self):
        self.assertEqual(unicycler.misc.strip_read_extensions('file.fasta'), 'file')
        self.assertEqual(unicycler.misc.strip_read_extensions('file.fasta.gz'), 'file')
        self.assertEqual(unicycler.misc.strip_read_extensions('path/to/file.fasta'), 'file')
        self.assertEqual(unicycler.misc.strip_read_extensions('path/to/file.fasta.gz'), 'file')

    def test_add_line_breaks_to_sequence(self):
        self.assertEqual(unicycler.misc.add_line_breaks_to_sequence('ATGCTGATGAAAATACC', 4),
                         'ATGC\nTGAT\nGAAA\nATAC\nC\n')
        self.assertEqual(unicycler.misc.add_line_breaks_to_sequence('ATGCTGATGAAAATACC', 8),
                         'ATGCTGAT\nGAAAATAC\nC\n')
        self.assertEqual(unicycler.misc.add_line_breaks_to_sequence('ATGCTGATGAAAATACC', 80),
                         'ATGCTGATGAAAATACC\n')

    def test_colour(self):
        t = 'test string'
        self.assertEqual(t, unicycler.misc.colour(t, 'normal'))
        self.assertEqual(t, unicycler.misc.colour(t, ''))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'green'))
        self.assertEqual(unicycler.misc.green(t), unicycler.misc.colour(t, 'green'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_green'))
        self.assertEqual(unicycler.misc.bold_green(t), unicycler.misc.colour(t, 'bold_green'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'red'))
        self.assertEqual(unicycler.misc.red(t), unicycler.misc.colour(t, 'red'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_red'))
        self.assertEqual(unicycler.misc.bold_red(t), unicycler.misc.colour(t, 'bold_red'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold'))
        self.assertEqual(unicycler.misc.bold(t), unicycler.misc.colour(t, 'bold'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_underline'))
        self.assertEqual(unicycler.misc.bold_underline(t),
                         unicycler.misc.colour(t, 'bold_underline'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'underline'))
        self.assertEqual(unicycler.misc.underline(t), unicycler.misc.colour(t, 'underline'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'dim'))
        self.assertEqual(unicycler.misc.dim(t), unicycler.misc.colour(t, 'dim'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'dim_underline'))
        self.assertEqual(unicycler.misc.dim_underline(t),
                         unicycler.misc.colour(t, 'dim_underline'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_yellow'))
        self.assertEqual(unicycler.misc.bold_yellow(t), unicycler.misc.colour(t, 'bold_yellow'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_yellow_underline'))
        self.assertEqual(unicycler.misc.bold_yellow_underline(t),
                         unicycler.misc.colour(t, 'bold_yellow_underline'))
        self.assertNotEqual(t, unicycler.misc.colour(t, 'bold_red_underline'))
        self.assertEqual(unicycler.misc.bold_red_underline(t),
                         unicycler.misc.colour(t, 'bold_red_underline'))

    def test_len_without_format(self):
        t = 'test string'
        self.assertEqual(len(t), unicycler.misc.len_without_format(unicycler.misc.green(t)))
        self.assertEqual(len(t), unicycler.misc.len_without_format(unicycler.misc.bold_red(t)))
        self.assertEqual(len(t),
                         unicycler.misc.len_without_format(unicycler.misc.bold_yellow_underline(t)))

    def test_remove_formatting(self):
        t = 'test string'
        self.assertEqual(t, unicycler.misc.remove_formatting(unicycler.misc.green(t)))
        self.assertEqual(t, unicycler.misc.remove_formatting(unicycler.misc.bold_red(t)))
        self.assertEqual(t,
                         unicycler.misc.remove_formatting(unicycler.misc.bold_yellow_underline(t)))

    def test_get_all_files_in_current_dir(self):
        starting_cwd = os.getcwd()
        test_dir = os.path.dirname(__file__)
        os.chdir(test_dir)
        self.assertTrue(os.path.basename(__file__) in unicycler.misc.get_all_files_in_current_dir())
        os.chdir(starting_cwd)

    def test_convert_fastq_to_fasta(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_1.fastq')
        test_fasta = os.path.join(os.path.dirname(__file__), 'temp_test.fasta')
        unicycler.misc.convert_fastq_to_fasta(test_fastq, test_fasta)
        fasta = unicycler.misc.load_fasta(test_fasta)
        self.assertEqual(len(fasta), 3)
        self.assertEqual(fasta[1][0], 'read_2')
        self.assertEqual(fasta[1][1], 'CACATACAGGCAGAGTGGCCGTGAAAGAAAGCAATCAGCGATGGTGCTCTGACGGGTTC'
                                      'GAGTTCTGCTGTGATAACGGAGAGAGACTGCGTGTCACGTTCGCGCTGGACTGCTGTGA'
                                      'TCGTGAG')
        os.remove(test_fasta)
