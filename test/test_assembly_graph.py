import unittest
import os
import unicycler.assembly_graph
import unicycler.misc


class TestAssemblyGraphFunctionsFastg(unittest.TestCase):
    """
    Tests various AssemblyGraph functions on a graph loaded from a SPAdes FASTG file.
    """

    def setUp(self):
        test_fastg = os.path.join(os.path.dirname(__file__), 'test.fastg')
        self.graph = unicycler.assembly_graph.AssemblyGraph(test_fastg, 25, paths_file=None,
                                                            insert_size_mean=401,
                                                            insert_size_deviation=60)

    def test_attributes(self):
        self.assertEqual(self.graph.overlap, 25)
        self.assertEqual(self.graph.insert_size_mean, 401)
        self.assertEqual(self.graph.insert_size_deviation, 60)

    def test_segments(self):
        self.assertEqual(len(self.graph.segments), 336)

    def test_forward_links(self):
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 904)

    def test_reverse_links(self):
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 904)

    def test_links_match(self):
        """
        Tests that each segment's forward links match up with its reverse complement's reverse
        links.
        """
        for segment, forward_link_list in self.graph.forward_links.items():
            reverse_link_list = self.graph.reverse_links[-segment]
            flipped_reverse_link_list = [-x for x in reverse_link_list]
            self.assertEqual(sorted(forward_link_list), sorted(flipped_reverse_link_list))

    def test_load_spades_paths(self):
        self.assertEqual(len(self.graph.paths), 0)
        paths_file = os.path.join(os.path.dirname(__file__), 'test.fastg.paths')
        self.graph.load_spades_paths(paths_file)
        self.assertEqual(len(self.graph.paths), 53)

    def test_get_median_read_depth(self):
        diff = abs(self.graph.get_median_read_depth() - 40.2)
        self.assertTrue(diff < 0.1)

    def test_normalise_read_depths(self):
        self.graph.normalise_read_depths()
        self.assertAlmostEqual(self.graph.get_median_read_depth(), 1.0)

    def test_get_total_length(self):
        self.assertEqual(self.graph.get_total_length(), 187896)

    def test_get_total_length_no_overlaps(self):
        self.assertEqual(self.graph.get_total_length_no_overlaps(), 179496)

    def test_total_dead_end_count(self):
        self.assertEqual(self.graph.total_dead_end_count(), 0)

    def test_dead_end_count(self):
        for seg_num in self.graph.segments:
            self.assertEqual(self.graph.dead_end_count(seg_num), 0)

    def test_save_to_fasta(self):
        temp_fasta = os.path.join(os.path.dirname(__file__), 'temp.fasta')
        self.graph.save_to_fasta(temp_fasta)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        self.assertEqual(fasta[0][0], '1')
        self.assertEqual(fasta[335][0], '336')
        self.assertEqual(fasta[0][1], '1 length=449 depth=82.60x')
        self.assertEqual(fasta[335][1], '336 length=185 depth=124.44x')
        self.assertTrue(fasta[0][2].startswith('ACCAGCCGCTGCGGGCCACCCGGAGCACGCGGCACATTGCCTTGATGCTG'
                                               'AACT'))
        self.assertTrue(fasta[335][2].endswith('TGTCGTGAAGCTTCACGAAGATGATTTTTTTGACGAAGAAGA'))
        os.remove(temp_fasta)

    def test_save_to_fasta_length(self):
        temp_fasta = os.path.join(os.path.dirname(__file__), 'temp.fasta')

        self.graph.save_to_fasta(temp_fasta, min_length=0)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=1)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=10)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=20)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=26)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 336)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=27)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 318)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=30)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 306)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=40)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 285)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=50)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 239)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=100)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 151)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=200)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 103)
        os.remove(temp_fasta)

        self.graph.save_to_fasta(temp_fasta, min_length=1000)
        fasta = unicycler.misc.load_fasta_with_full_header(temp_fasta)
        self.assertEqual(len(fasta), 40)
        os.remove(temp_fasta)

    def test_save_specific_segments_to_fasta(self):
        temp_fasta = os.path.join(os.path.dirname(__file__), 'temp.fasta')
        seg_nums = [1, 3, 5, 7, 9]
        segments = [self.graph.segments[x] for x in seg_nums]
        self.graph.save_specific_segments_to_fasta(temp_fasta, segments)
        fasta = unicycler.misc.load_fasta(temp_fasta)
        self.assertEqual(len(fasta), 5)
        self.assertEqual(fasta[0][0], '1')
        self.assertEqual(fasta[2][0], '5')
        self.assertTrue(fasta[0][1].startswith('ACCAGCCGCTGCGGGCCACCCGGAGCACGCGGCACATTGCCTTGATGCT'
                                               'GAACT'))
        self.assertTrue(fasta[2][1].endswith('TGCATGACAAAGTCATCGGGCATTATCTGAACATAAAACACTATCAATAAG'
                                             'TT'))
        os.remove(temp_fasta)

    def test_save_to_gfa(self):
        """
        This test saves the graph to GFA, loads the GFA graph and then compares the two graphs to
        make sure they are the same.
        """
        temp_gfa = os.path.join(os.path.dirname(__file__), 'temp.gfa')
        self.graph.save_to_gfa(temp_gfa, verbosity=0)
        graph2 = unicycler.assembly_graph.AssemblyGraph(temp_gfa, 25)
        self.assertEqual(self.graph.overlap, graph2.overlap)
        self.assertEqual(self.graph.insert_size_mean, graph2.insert_size_mean)
        self.assertEqual(self.graph.insert_size_deviation, graph2.insert_size_deviation)
        self.assertEqual(len(self.graph.segments), len(graph2.segments))
        link_count_1, link_count_2 = 0, 0
        for link_list in self.graph.forward_links.values():
            link_count_1 += len(link_list)
        for link_list in graph2.forward_links.values():
            link_count_2 += len(link_list)
        self.assertEqual(link_count_1, link_count_2)
        os.remove(temp_gfa)

    def test_get_all_gfa_link_lines(self):
        gfa_link_lines = self.graph.get_all_gfa_link_lines()
        self.assertEqual(gfa_link_lines.count('\n'), 452)
        self.assertEqual(gfa_link_lines.count('25M'), 452)

    def test_filter_by_read_depth(self):
        # A loop segment can be removed only when its depth drops below the threshold.
        self.assertEqual(len(self.graph.segments), 336)
        self.graph.filter_by_read_depth(0.5)
        self.assertEqual(len(self.graph.segments), 336)
        self.graph.segments[68].depth = 21.0
        self.graph.filter_by_read_depth(0.5)
        self.assertEqual(len(self.graph.segments), 336)
        self.graph.segments[68].depth = 20.0
        self.graph.filter_by_read_depth(0.5)
        self.assertEqual(len(self.graph.segments), 335)

        # A low-depth segment is only removed if its a dead end.
        self.graph.segments[306].depth = 0.1
        self.graph.filter_by_read_depth(0.5)
        self.assertEqual(len(self.graph.segments), 335)
        self.graph.remove_segments([273])
        self.assertEqual(len(self.graph.segments), 334)
        self.graph.filter_by_read_depth(0.5)
        self.assertEqual(len(self.graph.segments), 333)

    # def test_filter_homopolymer_loops(self):
    #     pass

    def test_remove_segments(self):
        self.assertEqual(len(self.graph.segments), 336)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 904)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 904)
        self.graph.remove_segments([276])
        self.assertEqual(len(self.graph.segments), 335)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 902)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 902)
        self.graph.remove_segments([273])
        self.assertEqual(len(self.graph.segments), 334)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 894)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 894)
        self.graph.remove_segments([67, 108, 222, 297])
        self.assertEqual(len(self.graph.segments), 330)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 870)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 870)

    def test_remove_small_components(self):
        self.assertEqual(len(self.graph.segments), 336)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 904)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 904)
        self.graph.remove_small_components(5000, 0)
        self.assertEqual(len(self.graph.segments), 336)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 904)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 904)
        self.graph.remove_small_components(6000, 0)
        self.assertEqual(len(self.graph.segments), 335)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 902)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 902)
        self.graph.remove_small_components(180000, 0)
        self.assertEqual(len(self.graph.segments), 335)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 902)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 902)
        self.graph.remove_small_components(190000, 0)
        self.assertEqual(len(self.graph.segments), 0)
        self.assertEqual(sum(len(x) for x in self.graph.forward_links.values()), 0)
        self.assertEqual(sum(len(x) for x in self.graph.reverse_links.values()), 0)

    # def test_remove_small_dead_ends(self):
    #     pass

    # def test_merge_all_possible(self):
    #     pass

    # def test_merge_simple_path(self):
    #     pass

    # def test_get_mean_path_depth(self):
    #     pass


# class TestAssemblyGraphFunctionsGfa(unittest.TestCase):
#     """
#     Tests various AssemblyGraph functions on a graph loaded from a GFA file.
#     """
#
#     def setUp(self):
#         pass
#
#     def test_segments(self):
#         pass
#
#     def test_forward_links(self):
#         pass
#
#     def test_reverse_links(self):
#         pass