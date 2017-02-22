"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This script tries to build Unicycler with various compilers and checks to see if they worked!

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import subprocess
import shutil
import os
import sys
import re

sys.path.insert(0, os.getcwd())
import unicycler.misc


def main():
    print()
    col_widths = [12, 27, 10, 23, 7, 7]
    header_row = ['Compiler', 'Location', 'Version', 'Build command', 'Build', 'Run']
    unicycler.misc.print_table([header_row], col_separation=2, header_format='underline', indent=0,
                               alignments='LLLLL', fixed_col_widths=col_widths)

    compiled_functions = os.path.join(os.getcwd(), 'unicycler', 'cpp_functions.so')

    for compiler in ['c++', 'g++', 'icpc', 'clang++',
                     'g++-4.8', 'g++-4.9', 'g++-5', 'g++-6',
                     'clang++-3.3', 'clang++-3.4', 'clang++-3.5', 'clang++-3.6', 'clang++-3.7',
                     'clang++-3.8', 'clang++-3.9']:
        clean_build()
        location = shutil.which(compiler)
        if location is None:
            location, version, make_command, build, run = 'not found', '', '' ,'', ''
        else:
            version = get_compiler_version(compiler)

            make_command = 'make CXX=' + compiler + ' -j'
            p = subprocess.Popen(make_command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, shell=True)
            p.communicate()

            if not os.path.isfile(compiled_functions):
                build = 'failed'
                run = ''
            else:
                build = 'success'

                # Make sure it actually works!
                if align_read():
                    run = 'success'
                else:
                    run = 'failed'

        unicycler.misc.print_table([[compiler, location, version, make_command, build, run]],
                                   col_separation=2, header_format='normal', indent=0,
                                   alignments='LLLLL', fixed_col_widths=col_widths,
                                   sub_colour={'success': 'green', 'failed': 'red',
                                               'not found': 'dim'})

    clean_build()


def clean_build():
    p = subprocess.Popen('make distclean', stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    p.communicate()


def align_read():
    ref_fasta = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fasta')
    read_fastq = os.path.join(os.path.dirname(__file__), 'test_semi_global_alignment.fastq')
    unicycler_align = os.path.join(os.getcwd(), 'unicycler_align-runner.py')

    temp_sam = 'TEMP_TEST_' + str(os.getpid()) + '.sam'

    p = subprocess.Popen([unicycler_align, '--ref', ref_fasta, '--reads', read_fastq,
                          '--sam', temp_sam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()

    # If all worked well, the SAM file should have 20 lines.
    if os.path.isfile(temp_sam):
        with open(temp_sam, 'rt') as sam:
            sam_line_count = sum(1 for _ in sam)
        os.remove(temp_sam)
    else:
        sam_line_count = 0
    return sam_line_count == 20


def get_compiler_version(compiler):
    p = subprocess.Popen(compiler + ' --version', stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, shell=True)
    stdout, _ = (x.decode().strip() for x in p.communicate())
    if 'version ' in stdout:
        version = stdout.split('version ')[-1].split()[0]
    else:
        try:
            version = re.findall('\d+\.\d+\.\d+', stdout)[0]
        except IndexError:
            version = '?'
    return version


if __name__ == '__main__':
    main()
