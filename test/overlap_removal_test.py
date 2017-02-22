"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This script generates random sequences with lots of repeats and then removes overlaps. It outputs
a table of information with the time taken for overlap removal.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import subprocess
import shutil
import random
import datetime
import sys

sys.path.insert(0, os.getcwd())
import unicycler.assembly_graph
import unicycler.misc
import unicycler.log
import test.fake_reads


def main():
    unicycler.log.logger = unicycler.log.Log(log_filename=None, stdout_verbosity_level=0)
    random.seed(0)
    print('\t'.join(['Sequence length', 'Repeat count', 'Segment count', 'Edge count',
                     'Overlap/kmer', 'Time (ms)']))
    try:
        while True:
            test_overlap_removal()

    # The user exits this script with Ctrl-C
    except KeyboardInterrupt:
        temp_dir = 'TEST_TEMP_' + str(os.getpid())
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir)


def test_overlap_removal():
    random_seq_length = random.randint(8, 20) ** 4
    repeat_count = random.randint(1, random_seq_length // 10)
    random_seq = make_repeaty_sequence(random_seq_length, repeat_count)
    out_dir = test.fake_reads.make_fake_reads(random_seq)

    output_line = [str(random_seq_length), str(repeat_count)]

    with open(os.path.join(out_dir, 'original_seq.fasta'), 'w') as fasta:
        fasta.write('>seq')
        fasta.write('\n')
        fasta.write(random_seq)
        fasta.write('\n')

    run_spades(out_dir)

    check_one_graph_overlap(out_dir, 21, output_line)
    check_one_graph_overlap(out_dir, 41, output_line)
    check_one_graph_overlap(out_dir, 61, output_line)
    check_one_graph_overlap(out_dir, 81, output_line)

    shutil.rmtree(out_dir)


def make_repeaty_sequence(length, repeat_count):
    seq = unicycler.misc.get_random_sequence(length)
    for i in range(repeat_count):
        repeat_length = random.randint(10, 500)
        repeat_instances = random.randint(2, 20)
        repeat_seq = unicycler.misc.get_random_sequence(repeat_length)
        for j in range(repeat_instances):
            if random.randint(0, 1) == 1:
                repeat_seq = unicycler.misc.reverse_complement(repeat_seq)
            repeat_pos = random.randint(0, length-repeat_length)
            seq = seq[:repeat_pos] + repeat_seq + seq[repeat_pos + repeat_length:]
            assert len(seq) == length
    return seq


def run_spades(out_dir):
    reads_1 = os.path.join(out_dir, 'reads_1.fastq')
    reads_2 = os.path.join(out_dir, 'reads_2.fastq')

    spades_cmd = ['spades.py', '-1', reads_1, '-2', reads_2, '-o', out_dir,
                  '--only-assembler', '-k', '21,41,61,81']
    p = subprocess.Popen(spades_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout.decode(), stderr.decode()


def check_one_graph_overlap(out_dir, k_size, output_line):
    k_dir = 'K' + str(k_size)
    graph_file = os.path.join(out_dir, k_dir, 'assembly_graph.fastg')
    graph = unicycler.assembly_graph.AssemblyGraph(graph_file, k_size)
    seg_count = len(graph.segments)
    edge_count = sum(len(x) for x in graph.forward_links.values()) // 2

    start_time = datetime.datetime.now()
    graph.remove_all_overlaps()
    end_time = datetime.datetime.now()
    milliseconds = (end_time - start_time).total_seconds() * 1000

    assert graph.overlap == 0
    graph.save_to_gfa(graph_file + '.gfa')
    print('\t'.join(output_line + [str(seg_count), str(edge_count), str(k_size),
                                   str(milliseconds)]))


if __name__ == '__main__':
    main()
