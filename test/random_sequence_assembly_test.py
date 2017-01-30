"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This script generates random sequences, shreds them, assembles them and makes sure that the output
sequences match the input sequences.

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
import unicycler.unicycler
import unicycler.assembly_graph
import unicycler.misc
import random
import datetime
import fake_reads


def main():
    random.seed(0)
    print('\t'.join(['Test type', 'Seq length', 'Time (ms)', 'Command']))
    while True:
        test_circular_no_repeat()
        test_circular_one_repeat()
        test_stdout_size()


def run_unicycler(out_dir, option_code, verbosity=None):
    """
    This function runs Unicycler. It uses different options, based on the iteration.
    """
    unicycler_runner = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    'unicycler-runner.py')
    reads_1 = os.path.join(out_dir, 'reads_1.fastq')
    reads_2 = os.path.join(out_dir, 'reads_2.fastq')

    unicycler_cmd = [unicycler_runner, '-1', reads_1, '-2', reads_2, '-o', out_dir]
    if option_code % 5 == 1:
        unicycler_cmd.append('--no_rotate')
    if option_code % 5 == 2:
        unicycler_cmd.append('--no_pilon')
    if option_code % 5 == 3:
        unicycler_cmd.append('--no_correct')
    if option_code % 5 == 4:
        unicycler_cmd += ['--no_rotate', '--no_pilon', '--no_correct']
    if verbosity is not None:
        unicycler_cmd += ['--verbosity', str(verbosity)]

    start_time = datetime.datetime.now()
    p = subprocess.Popen(unicycler_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    end_time = datetime.datetime.now()
    milliseconds = (end_time - start_time).total_seconds() * 1000

    cmd_string = ' '.join(unicycler_cmd)

    return stdout.decode(), stderr.decode(), cmd_string, milliseconds


def get_assembly_fasta_and_graph(out_dir):
    fasta = unicycler.misc.load_fasta(os.path.join(out_dir, 'assembly.fasta'))
    graph = unicycler.assembly_graph.AssemblyGraph(os.path.join(out_dir, 'assembly.gfa'),
                                                   0, paths_file=None)
    return fasta, graph


def sequence_matches_any_rotation(seq_1, seq_2):
    seq_1_rev_comp = unicycler.misc.reverse_complement(seq_1)
    for i in range(len(seq_1)):
        rotated_seq = seq_1[i:] + seq_1[:i]
        if seq_2 == rotated_seq:
            return True
        rotated_seq_rev_comp = seq_1_rev_comp[i:] + seq_1_rev_comp[:i]
        if seq_2 == rotated_seq_rev_comp:
            return True
    return False


def test_circular_no_repeat():
    random_seq_length = random.randint(8, 20) ** 4
    random_seq = unicycler.misc.get_random_sequence(random_seq_length)
    out_dir = fake_reads.make_fake_reads(random_seq)
    option_code = random.randint(0, 4)
    stdout, stderr, cmd_string, ms = run_unicycler(out_dir, option_code)
    assert bool(stderr) is False
    fasta, graph = get_assembly_fasta_and_graph(out_dir)
    assert len(fasta) == 1
    assembled_seq = fasta[0][1]
    assert sequence_matches_any_rotation(random_seq, assembled_seq)
    assert len(graph.segments) == 1
    assert len(graph.forward_links) == 2
    assert len(graph.reverse_links) == 2
    assert graph.segments[1].depth > 0.9
    assert graph.segments[1].depth < 1.1

    print('\t'.join(['circular_no_repeat', str(random_seq_length), str(ms), cmd_string]))
    shutil.rmtree(out_dir)


def test_circular_one_repeat():
    random_seq_length = random.randint(9, 20) ** 4
    repeat = unicycler.misc.get_random_sequence(500)
    non_repeat_length_1 = (random_seq_length - 1000) // 2
    non_repeat_length_2 = random_seq_length - 1000 - non_repeat_length_1
    seq_1 = unicycler.misc.get_random_sequence(non_repeat_length_1)
    seq_2 = unicycler.misc.get_random_sequence(non_repeat_length_2)
    random_seq = seq_1 + repeat + seq_2 + repeat
    out_dir = fake_reads.make_fake_reads(random_seq)
    option_code = random.randint(0, 4)
    stdout, stderr, cmd_string, ms = run_unicycler(out_dir, option_code)
    assert bool(stderr) is False
    fasta, graph = get_assembly_fasta_and_graph(out_dir)
    assert len(fasta) == 3
    seq_1 = fasta[0][1]
    seq_2 = fasta[1][1]
    seq_3 = fasta[2][1]

    assert len(graph.segments) == 3
    assert len(graph.forward_links) == 6
    assert len(graph.reverse_links) == 6
    assert graph.segments[1].depth > 0.9
    assert graph.segments[1].depth < 1.1
    assert graph.segments[2].depth > 0.9
    assert graph.segments[2].depth < 1.1
    assert graph.segments[3].depth > 1.9
    assert graph.segments[3].depth < 2.1

    assembled_len = len(seq_1) + len(seq_2) + 2 * len(seq_3)
    assert assembled_len == random_seq_length

    repeat_forward_links = set(graph.forward_links[3])
    if 1 in repeat_forward_links:
        s_1 = seq_1
    else:
        s_1 = unicycler.misc.reverse_complement(seq_1)
    if 2 in repeat_forward_links:
        s_2 = seq_2
    else:
        s_2 = unicycler.misc.reverse_complement(seq_2)
    assembled_seq = s_1 + seq_3 + s_2 + seq_3

    assert len(assembled_seq) == random_seq_length
    assert sequence_matches_any_rotation(random_seq, assembled_seq)

    print('\t'.join(['circular_one_repeat', str(random_seq_length), str(ms), cmd_string]))
    shutil.rmtree(out_dir)


def test_stdout_size():
    stdout_sizes = []
    for verbosity in range(4):
        random_seq_length = random.randint(1000, 5000)
        random_seq = unicycler.misc.get_random_sequence(random_seq_length)
        out_dir = fake_reads.make_fake_reads(random_seq)
        stdout, stderr, cmd_string, ms = run_unicycler(out_dir, 0, verbosity)
        assert bool(stderr) is False
        stdout_sizes.append(len(stdout))
        shutil.rmtree(out_dir)
        print('\t'.join(['stdout_size', str(random_seq_length), str(ms), cmd_string]))
    assert stdout_sizes[0] == 0
    assert stdout_sizes[1] > 0
    assert stdout_sizes[2] >= stdout_sizes[1]
    assert stdout_sizes[3] >= stdout_sizes[2]


if __name__ == '__main__':
    main()
