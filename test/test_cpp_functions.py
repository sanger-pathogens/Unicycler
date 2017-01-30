import unittest
import os
import unicycler.cpp_function_wrappers
import unicycler.read_ref
import unicycler.alignment

class TestSemiGlobalAlignment(unittest.TestCase):
    pass


class TestFullyGlobalAlignment(unittest.TestCase):
    pass


class TestPathAlignment(unittest.TestCase):
    pass


class TestMultipleSequenceAlignment(unittest.TestCase):

    def setUp(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_2.fastq')
        self.reads, _, _ = unicycler.read_ref.load_long_reads(test_fastq, 0)
        self.scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        read_count = len(self.reads)
        self.seqs = [self.reads[str(i+1)].sequence for i in range(read_count)]
        self.quals = [self.reads[str(i+1)].qualities for i in range(read_count)]

    def test_consensus_1(self):
        seqs = self.seqs[:3]
        quals = self.quals[:3]

        print(seqs)  # TEMP
        print(quals)  # TEMP

        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, 'ACGTGTGTTCGCCATGGCTGCGATTTCT')
