"""
This script generates random sequences with lots of repeats and then removes overlaps.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
import shutil
import unicycler.assembly_graph
import unicycler.misc
import random
import datetime


def main():
    random.seed(0)
    print('\t'.join(['Sequence length', 'Repeat count', 'Segment count', 'Edge count',
                     'Overlap/kmer', 'Time (ms)']))
    while True:
        test_overlap_removal()


def test_overlap_removal():
    random_seq_length = random.randint(5000, 50000)
    repeat_count = random.randint(1, random_seq_length // 10)
    random_seq = make_repeaty_sequence(random_seq_length, repeat_count)
    out_dir = make_fake_reads(random_seq)

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
    graph.remove_all_overlaps(0)
    end_time = datetime.datetime.now()
    milliseconds = (end_time - start_time).total_seconds() * 1000

    assert graph.overlap == 0
    graph.save_to_gfa(graph_file + '.gfa', 0)
    print('\t'.join(output_line + [str(seg_count), str(edge_count), str(k_size),
                                   str(milliseconds)]))


def make_fake_qual_string(length):
    qual_string = ''
    for i in range(length):
        qual_string += chr(random.randint(70, 75))
    return qual_string


def make_fake_reads(seq):
    """
    This function makes super-simple fake reads (no errors and even distribution) for a circular
    genome and saves them to FASTQ files.
    """
    read_length = 100
    insert_length = 300
    looped_seq = seq + seq[0:insert_length]

    out_dir = 'TEST_TEMP_' + str(os.getpid())
    os.makedirs(out_dir)

    reads_1 = os.path.join(out_dir, 'reads_1.fastq')
    reads_2 = os.path.join(out_dir, 'reads_2.fastq')

    read_num = 1
    with open(reads_1, 'wt') as r_1, open(reads_2, 'wt') as r_2:
        for i in range(0, 2):
            for j in range(0, len(looped_seq) - insert_length + 1):
                if i == 1:
                    looped_seq = unicycler.misc.reverse_complement(looped_seq)
                insert = looped_seq[j:j+insert_length]
                read_1 = insert[:read_length]
                read_2 = insert[-read_length:]
                r_1.write('@read_' + str(read_num) + '/1\n')
                r_2.write('@read_' + str(read_num) + '/2\n')
                r_1.write(read_1 + '\n')
                r_2.write(read_2 + '\n')
                r_1.write('+\n')
                r_2.write('+\n')
                r_1.write(make_fake_qual_string(read_length) + '\n')
                r_2.write(make_fake_qual_string(read_length) + '\n')
                read_num += 1
    return out_dir


if __name__ == '__main__':
    main()