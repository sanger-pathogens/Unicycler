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
import subprocess
import shutil
import re
import stat


class TestDependencies(unittest.TestCase):
    def setUp(self):
        sample_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sample_data')
        self.reads_1 = os.path.join(sample_dir, 'short_reads_1.fastq.gz')
        self.reads_2 = os.path.join(sample_dir, 'short_reads_2.fastq.gz')
        self.reads_long = os.path.join(sample_dir, 'long_reads_low_depth.fastq.gz')
        self.out_dir = 'TEMP_' + str(os.getpid())

    def tearDown(self):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

    def run_unicycler(self, options):
        unicycler_runner = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                        'unicycler-runner.py')
        unicycler_cmd = [unicycler_runner, '-1', self.reads_1, '-2', self.reads_2,
                         '-l', self.reads_long, '-o', self.out_dir] + options
        p = subprocess.Popen(unicycler_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        return stdout.decode(), stderr.decode(), p.returncode

    def test_spades_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--spades_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'spades.py\s+not found', stdout)))
        self.assertTrue('could not find SPAdes' in stderr)
        self.assertEqual(ret_code, 1)

    def test_makeblastdb_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--makeblastdb_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'makeblastdb\s+not found', stdout)))
        self.assertTrue('could not find makeblastdb' in stderr)
        self.assertEqual(ret_code, 1)

    def test_tblastn_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--tblastn_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'tblastn\s+not found', stdout)))
        self.assertTrue('could not find tblastn' in stderr)
        self.assertEqual(ret_code, 1)

    def test_bowtie2_build_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--bowtie2_build_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'bowtie2-build\s+not found', stdout)))
        self.assertTrue('could not find bowtie2-build' in stderr)
        self.assertEqual(ret_code, 1)

    def test_bowtie2_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--bowtie2_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'bowtie2\s+not found', stdout)))
        self.assertTrue('could not find bowtie2' in stderr)
        self.assertEqual(ret_code, 1)

    def test_samtools_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--samtools_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'samtools\s+not found', stdout)))
        self.assertTrue('could not find samtools' in stderr)
        self.assertEqual(ret_code, 1)

    def test_java_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--java_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'java\s+not found', stdout)))
        self.assertTrue('could not find java' in stderr)
        self.assertEqual(ret_code, 1)

    def test_pilon_not_found(self):
        stdout, stderr, ret_code = self.run_unicycler(['--pilon_path', 'not_a_real_path'])
        self.assertTrue(bool(re.search(r'pilon\s+not found', stdout)))
        self.assertTrue('could not find pilon')
        self.assertEqual(ret_code, 1)

    def test_bad_pilon(self):
        bad_pilon = 'TEMP_bad_pilon_' + str(os.getpid())
        open(bad_pilon, 'a').close()
        os.chmod(bad_pilon, os.stat(bad_pilon).st_mode | stat.S_IEXEC)
        stdout, stderr, ret_code = self.run_unicycler(['--pilon_path', bad_pilon])
        self.assertTrue(bool(re.search(r'pilon\s+bad', stdout)))
        self.assertTrue('Pilon was found' in stderr)
        self.assertTrue('but does not work' in stderr)
        self.assertEqual(ret_code, 1)
        os.remove(bad_pilon)

    def test_no_rotate(self):
        stdout, stderr, ret_code = self.run_unicycler(['--spades_path', 'not_a_real_path',
                                                       '--no_rotate'])
        self.assertTrue(bool(re.search(r'spades.py\s+not found', stdout)))
        self.assertTrue(bool(re.search(r'makeblastdb\s+not used', stdout)))
        self.assertTrue(bool(re.search(r'tblastn\s+not used', stdout)))

    def test_no_pilon(self):
        stdout, stderr, ret_code = self.run_unicycler(['--spades_path', 'not_a_real_path',
                                                       '--no_pilon'])
        self.assertTrue(bool(re.search(r'spades.py\s+not found', stdout)))
        self.assertTrue(bool(re.search(r'bowtie2-build\s+not used', stdout)))
        self.assertTrue(bool(re.search(r'bowtie2\s+not used', stdout)))
        self.assertTrue(bool(re.search(r'samtools\s+not used', stdout)))
        self.assertTrue(bool(re.search(r'java\s+not used', stdout)))
        self.assertTrue(bool(re.search(r'pilon\s+not used', stdout)))

    def test_verbosity_1(self):
        stdout, stderr, ret_code = self.run_unicycler(['--spades_path', 'not_a_real_path',
                                                       '--verbosity', '1'])
        self.assertTrue(bool(re.search(r'spades.py\s+not found', stdout)))
        self.assertTrue(bool(re.search(r'Program\s+Version\s+Status', stdout)))
        self.assertFalse(bool(re.search(r'Program\s+Version\s+Status\s+Path', stdout)))

    def test_verbosity_2(self):
        stdout, stderr, ret_code = self.run_unicycler(['--spades_path', 'not_a_real_path',
                                                       '--verbosity', '2'])
        self.assertTrue(bool(re.search(r'spades.py\s+not found', stdout)))
        self.assertTrue(bool(re.search(r'Program\s+Version\s+Status\s+Path', stdout)))
