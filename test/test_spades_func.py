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
import unicycler.spades_func


class TestSPAdesFunc(unittest.TestCase):

    def test_get_read_lengths_1(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_misc.fastq')
        read_lengths = sorted(unicycler.spades_func.get_read_lengths(test_fastq))
        self.assertEqual(read_lengths, [125, 125, 125])

    def test_get_read_lengths_2(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fastq')
        read_lengths = sorted(unicycler.spades_func.get_read_lengths(test_fastq))
        self.assertEqual(read_lengths, [100, 150, 200, 300, 300, 300, 300, 300, 300])

    def test_get_read_count_1(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_misc.fastq')
        self.assertEqual(unicycler.spades_func.get_read_count(test_fastq), 3)

    def test_get_read_count_2(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fastq')
        self.assertEqual(unicycler.spades_func.get_read_count(test_fastq), 9)

    def test_bad_fastq_1(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_bad_reads_1.fastq')
        with self.assertRaises(unicycler.spades_func.BadFastq):
            unicycler.spades_func.get_read_count(test_fastq)

    def test_bad_fastq_2(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_bad_reads_2.fastq')
        with self.assertRaises(unicycler.spades_func.BadFastq):
            unicycler.spades_func.get_read_count(test_fastq)

    def test_build_spades_command_1(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             0, '1.fq.gz', '2.fq.gz', None, True,
                                                             False, None)
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21', '--threads', '16',
                                   '--isolate', '-1', '1.fq.gz', '-2', '2.fq.gz', '-m', '1024'])

    def test_build_spades_command_2(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             1, '1.fq.gz', '2.fq.gz', None, True,
                                                             False, None)
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21,31', '--threads', '16',
                                   '--restart-from', 'k21', '-m', '1024'])

    def test_build_spades_command_3(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             2, '1.fq.gz', '2.fq.gz', None, True,
                                                             False, None)
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21,31,41', '--threads', '16',
                                   '--restart-from', 'k31', '-m', '1024'])

    def test_build_spades_command_4(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             0, None, None, 's.fq.gz', False, True,
                                                             None)
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21', '--threads', '16',
                                   '--isolate', '-s', 's.fq.gz', '-m', '1024'])

    def test_build_spades_command_5(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             0, '1.fq.gz', '2.fq.gz', None, True,
                                                             False, '-m 123')
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21', '--threads', '16',
                                   '--isolate', '-1', '1.fq.gz', '-2', '2.fq.gz', '-m', '123'])

    def test_build_spades_command_6(self):
        command = unicycler.spades_func.build_spades_command('spades.py', 'out', 16, [21, 31, 41],
                                                             0, '1.fq.gz', '2.fq.gz', None, True,
                                                             False, '--tmp-dir abc')
        self.assertEqual(command, ['spades.py', '-o', 'out', '-k', '21', '--threads', '16',
                                   '--isolate', '-1', '1.fq.gz', '-2', '2.fq.gz', '--tmp-dir',
                                   'abc', '-m', '1024'])

