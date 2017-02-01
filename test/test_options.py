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
import unicycler.unicycler
import sys
import os
import io
from contextlib import contextmanager


@contextmanager
def captured_output():
    new_out, new_err = io.StringIO(), io.StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


def get_default_from_help(option, help_text):
    return help_text.split(option)[1].split('default:')[1].split(')')[0].strip()


class TestUnicyclerOptions(unittest.TestCase):

    def test_no_options(self):
        sys.argv = [sys.argv[0]]
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit):
                unicycler.unicycler.get_arguments()
        stdout = out.getvalue().strip()
        stderr = err.getvalue().strip()
        self.assertEqual(len(stdout), 0)
        self.assertTrue('error' in stderr)
        self.assertTrue('required' in stderr)

    def test_absolute_paths(self):
        reads_1, reads_2 = '/path/to/reads_1.fastq.gz', '/path/to/reads_1.fastq.gz'
        unpaired, long = '/path/to/unpaired.fastq.gz', '/path/to/long.fastq.gz'
        out_dir = '/path/to/out_dir'
        sys.argv = [sys.argv[0]] + \
                   ['-1', reads_1, '-2', reads_2, '-s', unpaired, '-l', long, '-o', out_dir]
        with captured_output() as (out, err):
            args = unicycler.unicycler.get_arguments()
        stdout = out.getvalue().strip()
        stderr = err.getvalue().strip()
        self.assertEqual(len(stdout), 0)
        self.assertEqual(len(stderr), 0)
        self.assertEqual(args.short1, reads_1)
        self.assertEqual(args.short2, reads_2)
        self.assertEqual(args.unpaired, unpaired)
        self.assertEqual(args.long, long)
        self.assertEqual(args.out, out_dir)

    def test_relative_paths(self):
        reads_1, reads_2 = 'reads_1.fastq.gz', 'reads_1.fastq.gz'
        unpaired, long = 'unpaired.fastq.gz', 'long.fastq.gz'
        out_dir = 'out_dir'
        sys.argv = [sys.argv[0]] + \
                   ['-1', reads_1, '-2', reads_2, '-s', unpaired, '-l', long, '-o', out_dir]
        with captured_output() as (out, err):
            args = unicycler.unicycler.get_arguments()
        stdout = out.getvalue().strip()
        stderr = err.getvalue().strip()
        self.assertEqual(len(stdout), 0)
        self.assertEqual(len(stderr), 0)
        reads_1_full_path = os.path.join(os.getcwd(), reads_1)
        reads_2_full_path = os.path.join(os.getcwd(), reads_2)
        unpaired_full_path = os.path.join(os.getcwd(), unpaired)
        long_full_path = os.path.join(os.getcwd(), long)
        out_dir_full_path = os.path.join(os.getcwd(), out_dir)
        self.assertEqual(args.short1, reads_1_full_path)
        self.assertEqual(args.short2, reads_2_full_path)
        self.assertEqual(args.unpaired, unpaired_full_path)
        self.assertEqual(args.long, long_full_path)
        self.assertEqual(args.out, out_dir_full_path)

    def test_defaults(self):
        sys.argv = [sys.argv[0], '--help_all']
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit):
                unicycler.unicycler.get_arguments()
        help_text = out.getvalue().strip().split('Unicycler: a hybrid assembly pipeline for '
                                                 'bacterial genomes')[1]
        verbosity_default = int(get_default_from_help('--verbosity', help_text))
        min_fasta_length_default = int(get_default_from_help('--min_fasta_length', help_text))
        keep_temp_default = int(get_default_from_help('--keep_temp', help_text))
        threads_default = int(get_default_from_help('--threads', help_text))
        # mode_default = get_default_from_help('--mode', help)
        expected_linear_seqs_default = int(get_default_from_help('--expected_linear_seqs',
                                                                 help_text))
        spades_path_default = get_default_from_help('--spades_path', help_text)
        min_kmer_frac_default = float(get_default_from_help('--min_kmer_frac', help_text))
        max_kmer_frac_default = float(get_default_from_help('--max_kmer_frac', help_text))
        kmer_count_default = int(get_default_from_help('--kmer_count', help_text))
        start_gene_id_default = float(get_default_from_help('--start_gene_id', help_text))
        start_gene_cov_default = float(get_default_from_help('--start_gene_cov', help_text))
        makeblastdb_path_default = get_default_from_help('--makeblastdb_path', help_text)
        tblastn_path_default = get_default_from_help('--tblastn_path', help_text)
        bowtie2_path_default = get_default_from_help('--bowtie2_path', help_text)
        bowtie2_build_path_default = get_default_from_help('--bowtie2_build_path', help_text)
        samtools_path_default = get_default_from_help('--samtools_path', help_text)
        pilon_path_default = get_default_from_help('--pilon_path', help_text)
        java_path_default = get_default_from_help('--java_path', help_text)
        min_polish_size_default = int(get_default_from_help('--min_polish_size', help_text))
        min_component_size_default = int(get_default_from_help('--min_component_size', help_text))
        min_dead_end_size_default = int(get_default_from_help('--min_dead_end_size', help_text))
        scores_default = get_default_from_help('--scores', help_text)

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.verbosity, verbosity_default)
        self.assertEqual(args.min_fasta_length, min_fasta_length_default)
        self.assertEqual(args.keep_temp, keep_temp_default)
        self.assertEqual(args.threads, threads_default)
        self.assertEqual(args.expected_linear_seqs, expected_linear_seqs_default)
        self.assertEqual(args.spades_path, spades_path_default)
        self.assertEqual(args.min_kmer_frac, min_kmer_frac_default)
        self.assertEqual(args.max_kmer_frac, max_kmer_frac_default)
        self.assertEqual(args.kmer_count, kmer_count_default)
        self.assertEqual(args.start_gene_id, start_gene_id_default)
        self.assertEqual(args.start_gene_cov, start_gene_cov_default)
        self.assertEqual(args.makeblastdb_path, makeblastdb_path_default)
        self.assertEqual(args.tblastn_path, tblastn_path_default)
        self.assertEqual(args.bowtie2_path, bowtie2_path_default)
        self.assertEqual(args.bowtie2_build_path, bowtie2_build_path_default)
        self.assertEqual(args.samtools_path, samtools_path_default)
        self.assertEqual(args.pilon_path, pilon_path_default)
        self.assertEqual(args.java_path, java_path_default)
        self.assertEqual(args.min_polish_size, min_polish_size_default)
        self.assertEqual(args.min_component_size, min_component_size_default)
        self.assertEqual(args.min_dead_end_size, min_dead_end_size_default)
        self.assertEqual(args.scores, scores_default)

    def test_modes(self):
        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'conservative']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 0)
        conservative_min_bridge_qual = args.min_bridge_qual

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'normal']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 1)
        normal_min_bridge_qual = args.min_bridge_qual

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'bold']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 2)
        bold_min_bridge_qual = args.min_bridge_qual

        self.assertTrue(conservative_min_bridge_qual > normal_min_bridge_qual)
        self.assertTrue(normal_min_bridge_qual > bold_min_bridge_qual)

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'conservative', '--min_bridge_qual', '98.7']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 0)
        self.assertEqual(args.min_bridge_qual, 98.7)

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'normal', '--min_bridge_qual', '54.3']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 1)
        self.assertEqual(args.min_bridge_qual, 54.3)

        sys.argv = [sys.argv[0], '-1', 'reads_1.fastq.gz', '-2', 'reads_2.fastq.gz',
                    '-o', 'out_dir', '--mode', 'bold', '--min_bridge_qual', '12.3']
        with captured_output() as (_, _):
            args = unicycler.unicycler.get_arguments()
        self.assertEqual(args.mode, 2)
        self.assertEqual(args.min_bridge_qual, 12.3)
