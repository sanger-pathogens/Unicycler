import unittest
import os
import subprocess
import shutil
import unicycler.unicycler
import unicycler.assembly_graph
import unicycler.misc
import random
import test.fake_reads


REPLICATE_COUNT = 5


def run_unicycler(out_dir, option_code, verbosity=None):
    """
    This function runs Unicycler. It uses different options, based on the iteration.
    """
    unicycler_runner = os.path.join(os.path.dirname(os.path.dirname(__file__)),
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

    p = subprocess.Popen(unicycler_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout.decode(), stderr.decode()


def get_assembly_fasta_and_graph(out_dir):
    fasta = unicycler.misc.load_fasta(os.path.join(out_dir, 'assembly.fasta'))
    graph = unicycler.assembly_graph.AssemblyGraph(os.path.join(out_dir, 'assembly.gfa'),
                                                   0, paths_file=None)
    return fasta, graph


def sequence_matches_either_strand(seq_1, seq_2):
    if seq_1 == seq_2:
        return True
    else:
        return unicycler.misc.reverse_complement(seq_1) == seq_2


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


class TestRandomSeqAssemblies(unittest.TestCase):
    """
    Tests various AssemblyGraph functions on a graph loaded from a SPAdes FASTG file.
    """
    def test_1000_bp_circular_no_repeat(self):
        for i in range(REPLICATE_COUNT):
            random.seed(i)
            random_seq = unicycler.misc.get_random_sequence(1000)
            out_dir = test.fake_reads.make_fake_reads(random_seq)
            stdout, stderr = run_unicycler(out_dir, i)
            self.assertFalse(bool(stderr), msg=stderr)
            fasta, graph = get_assembly_fasta_and_graph(out_dir)
            self.assertEqual(len(fasta), 1)
            assembled_seq = fasta[0][1]
            self.assertTrue(sequence_matches_any_rotation(random_seq, assembled_seq))
            self.assertEqual(len(graph.segments), 1)
            self.assertEqual(len(graph.forward_links), 2)
            self.assertEqual(len(graph.reverse_links), 2)
            self.assertTrue(graph.segments[1].depth > 0.9)
            self.assertTrue(graph.segments[1].depth < 1.1)
            shutil.rmtree(out_dir)

    def test_5000_bp_circular_no_repeat(self):
        for i in range(REPLICATE_COUNT):
            random.seed(i)
            random_seq = unicycler.misc.get_random_sequence(5000)
            out_dir = test.fake_reads.make_fake_reads(random_seq)
            stdout, stderr = run_unicycler(out_dir, i)
            self.assertFalse(bool(stderr), msg=stderr)
            fasta, graph = get_assembly_fasta_and_graph(out_dir)
            self.assertEqual(len(fasta), 1)
            assembled_seq = fasta[0][1]
            self.assertTrue(sequence_matches_any_rotation(random_seq, assembled_seq))
            self.assertEqual(len(graph.segments), 1)
            self.assertEqual(len(graph.forward_links), 2)
            self.assertEqual(len(graph.reverse_links), 2)
            self.assertTrue(graph.segments[1].depth > 0.9)
            self.assertTrue(graph.segments[1].depth < 1.1)
            shutil.rmtree(out_dir)

    def test_5000_bp_circular_one_repeat(self):
        for i in range(REPLICATE_COUNT):
            random.seed(i)
            repeat = unicycler.misc.get_random_sequence(500)
            seq_1 = unicycler.misc.get_random_sequence(2500)
            seq_2 = unicycler.misc.get_random_sequence(1500)
            random_seq = seq_1 + repeat + seq_2 + repeat
            out_dir = test.fake_reads.make_fake_reads(random_seq)
            stdout, stderr = run_unicycler(out_dir, i)
            self.assertFalse(bool(stderr), msg=stderr)
            fasta, graph = get_assembly_fasta_and_graph(out_dir)
            self.assertEqual(len(fasta), 3)
            seq_1 = fasta[0][1]
            seq_2 = fasta[1][1]
            seq_3 = fasta[2][1]

            self.assertEqual(len(graph.segments), 3)
            self.assertEqual(len(graph.forward_links), 6)
            self.assertEqual(len(graph.reverse_links), 6)
            self.assertTrue(graph.segments[1].depth > 0.9)
            self.assertTrue(graph.segments[1].depth < 1.1)
            self.assertTrue(graph.segments[2].depth > 0.9)
            self.assertTrue(graph.segments[2].depth < 1.1)
            self.assertTrue(graph.segments[3].depth > 1.9)
            self.assertTrue(graph.segments[3].depth < 2.1)

            assembled_len = len(seq_1) + len(seq_2) + 2 * len(seq_3)
            self.assertEqual(assembled_len, 5000)

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

            self.assertEqual(len(assembled_seq), 5000)
            self.assertTrue(sequence_matches_any_rotation(random_seq, assembled_seq))

            shutil.rmtree(out_dir)

    def test_stdout_size(self):
        stdout_sizes = []
        for verbosity in range(4):
            random.seed(0)
            random_seq = unicycler.misc.get_random_sequence(1000)
            out_dir = test.fake_reads.make_fake_reads(random_seq)
            stdout, stderr = run_unicycler(out_dir, 0, verbosity)
            self.assertFalse(bool(stderr), msg=stderr)
            stdout_sizes.append(len(stdout))
            shutil.rmtree(out_dir)
        self.assertTrue(stdout_sizes[0] == 0)
        self.assertTrue(stdout_sizes[1] > 0)
        self.assertTrue(stdout_sizes[2] >= stdout_sizes[1])
        self.assertTrue(stdout_sizes[3] >= stdout_sizes[2])
