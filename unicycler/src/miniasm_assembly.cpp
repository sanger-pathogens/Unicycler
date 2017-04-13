// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "miniasm_assembly.h"

#include <iostream>
#include <fstream>

//#include <unistd.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <vector>
//#include <sys/stat.h>
#include "miniasm/sys.h"
#include "miniasm/sdict.h"
#include "miniasm/miniasm.h"
//#include "tclap/CmdLine.h"
//
//
//#define VERSION "0.3-unicycler"
//
//using namespace TCLAP;
using namespace std;


void miniasmAssembly(char * reads, char * overlaps, char * outputDir, int min_dp) {
    string paf_filename(overlaps);     // Input PAF mapping
    string reads_filename(reads);      // Input long reads
    string outdir(outputDir);          // Output directory

    bool prefilter_contained = false;  // prefilter clearly contained reads (2-pass required)
    int min_match = 100;               // min match length
    int min_span = 2000;               // min span
    int min_ovlp = 2000;               //
    int max_hang = 1000;               // max over hang length
    int gap_fuzz = 1000;               // max gap differences between reads for trans-reduction
    int bub_dist = 50000;              // max distance for bubble popping
    int max_ext = 4;                   // small unitig threshold
    int n_rounds = 2;                  // rounds of short overlap removal
    float min_iden = 0.05;             // min identity
    float int_frac = 0.8;              // min end-to-end match ratio
    float min_drop = 0.5;              // min overlap drop ratio
    float max_drop = 0.7;              // max overlap drop ratio
    float final_drop = 0.8;            // aggressive overlap drop ratio in the end
    bool bi_dir = false;               // both directions of an arc are present in input

    float cov = 40.0;
    string outfmt = "ug";

    string raw_string_graph = outdir + "/03_raw_string_graph.gfa";
    string transitive_reduction_string_graph = outdir + "/04_after_transitive_reduction.gfa";
    string tip_cut_string_graph = outdir + "/05_after_tip_cutting.gfa";
    string bubble_pop_string_graph = outdir + "/06_after_bubble_popping.gfa";
    string cut_overlaps_string_graph_1 = outdir + "/07_after_cutting_short_overlaps.gfa";
    string remove_internal_string_graph = outdir + "/08_after_removing_internal_seqs.gfa";
    string cut_overlaps_string_graph_2 = outdir + "/09_after_cutting_short_overlaps.gfa";
    string final_string_graph = outdir + "/10_final_string_graph.gfa";
    string miniasm_output = outdir + "/miniasm.out";
    string contained_read_list = outdir + "/contained_reads.txt";
    string all_read_list = outdir + "/all_reads.txt";

    // Redirect miniasm's output to a stringstream, instead of outputting it to stderr.
    // http://stackoverflow.com/questions/5419356/redirect-stdout-stderr-to-a-string
    std::ofstream outFile;
    outFile.open(miniasm_output);
    std::streambuf * old = std::cerr.rdbuf(outFile.rdbuf());

    sys_init();

    // If the user ran miniasm with -R, this extra step is carried out to remove contained reads.
    sdict_t *excluded_reads = 0;
    if (prefilter_contained) {
        cerr << "===> Step 0: removing contained reads <===\n";
        excluded_reads = prefilter_contained_reads(paf_filename.c_str(), min_span, min_match, max_hang, int_frac);
    }

    // Load in the PAF file of read-to-read alignments, excluding alignments which are too short.
    // If step 0 was run, contained reads are also excluded. Also a single PAF line is loaded in
    // both directions (unless the -b option was used which implies that both directions are
    // present in the PAF file).
    cerr << "===> Step 1: reading read mappings <===\n";
    size_t num_hits;
    sdict_t *read_dict = init_seq_dict();
    ma_hit_t * hits = read_hits_file(paf_filename.c_str(), min_span, min_match, read_dict, &num_hits, !bi_dir, excluded_reads);
    std::cerr << "\n";

    ma_sub_t *subreads = 0;

    cerr << "===> Step 2: 1-pass (crude) read selection <===\n";

    // Toss out reads which fail to meet the read depth threshold. It creates a sub object
    // which stores the start and end positions of a read which have met the min depth. This
    // first pass looks at the entire mappings (not clipped off at all).
    subreads = filter_reads_using_depth(min_dp, min_iden, 0, num_hits, hits, read_dict);

    // Toss out hits which fail to meet the minimum span threshold.
    num_hits = filter_hits_using_span(subreads, min_span, num_hits, hits);

    // Toss out hits which have too much overhang.
    num_hits = filter_hits_using_overhang(subreads, int(max_hang * 1.5), int(min_ovlp * 0.5), num_hits, hits, &cov);
    std::cerr << "\n";

    cerr << "===> Step 3: 2-pass (fine) read selection <===\n";
    ma_sub_t *subreads_2;

    // Toss out reads which fail to meet the read depth threshold again.
    // This second pass trims the mappings down making it harder to meet the threshold (i.e.
    // the trimmed mapping must meet the depth threshold).
    subreads_2 = filter_reads_using_depth(min_dp, min_iden, min_span/2, num_hits, hits, read_dict);

    // Filter hits using the minimum span threshold again. Since we just did a more stringent
    // run through filter_reads_using_depth, this can toss out more alignments than our first call to this
    // function.
    num_hits = filter_hits_using_span(subreads_2, min_span, num_hits, hits);

    merge_subreads(read_dict->n_seq, subreads, subreads_2);
    free(subreads_2);
    save_read_names(num_hits, hits, read_dict, subreads, all_read_list);

    // Toss out contained reads (this is a big one and gets rid of a lot).
    num_hits = remove_contained_reads(max_hang, int_frac, min_ovlp, read_dict, subreads, num_hits, hits, contained_read_list);
    std::cerr << "\n";

    hits = (ma_hit_t*)realloc(hits, num_hits * sizeof(ma_hit_t));
    asg_t * string_graph = 0;

    cerr << "===> Step 4: graph cleaning <===\n";
    string_graph = make_string_graph(max_hang, int_frac, min_ovlp, read_dict, subreads, num_hits, hits);
    save_string_graph(string_graph, read_dict, subreads, raw_string_graph, reads_filename.c_str());
    std::cerr << "\n";

    cerr << "===> Step 4.1: transitive reduction <===\n";
    asg_arc_del_trans(string_graph, gap_fuzz);
    save_string_graph(string_graph, read_dict, subreads, transitive_reduction_string_graph, reads_filename.c_str());
    std::cerr << "\n";

    cerr << "===> Step 4.2: initial tip cutting and bubble popping <===\n";
    cut_tips(string_graph, max_ext);
    save_string_graph(string_graph, read_dict, subreads, tip_cut_string_graph, reads_filename.c_str());
    pop_bubbles(string_graph, bub_dist);
    save_string_graph(string_graph, read_dict, subreads, bubble_pop_string_graph, reads_filename.c_str());
    std::cerr << "\n";

    cerr << "===> Step 4.3: cutting short overlaps (%d rounds in total) <===\n";
    for (int i = 0; i <= n_rounds; ++i) {
        float r = min_drop + (max_drop - min_drop) / n_rounds * i;
        if (asg_arc_del_short(string_graph, r) != 0) {
            cut_tips(string_graph, max_ext);
            pop_bubbles(string_graph, bub_dist);
        }
    }
    save_string_graph(string_graph, read_dict, subreads, cut_overlaps_string_graph_1, reads_filename.c_str());
    std::cerr << "\n";

    cerr << "===> Step 4.4: removing short internal sequences and bi-loops <===\n";
    cut_short_internal(string_graph, 1);
    cut_biloops(string_graph, max_ext);
    cut_tips(string_graph, max_ext);
    pop_bubbles(string_graph, bub_dist);
    save_string_graph(string_graph, read_dict, subreads, remove_internal_string_graph, reads_filename.c_str());
    std::cerr << "\n";

    cerr << "===> Step 4.5: aggressively cutting short overlaps <===\n";
    if (asg_arc_del_short(string_graph, final_drop) != 0) {
        cut_tips(string_graph, max_ext);
        pop_bubbles(string_graph, bub_dist);
    }
    save_string_graph(string_graph, read_dict, subreads, cut_overlaps_string_graph_2, reads_filename.c_str());
    std::cerr << "\n";

    save_string_graph(string_graph, read_dict, subreads, final_string_graph, reads_filename.c_str());
    destroy_string_graph(string_graph);

    // Clean up!
    free(subreads); free(hits);
    destroy_seq_dict(read_dict);
    if (excluded_reads)
        destroy_seq_dict(excluded_reads);

    cerr << "Real time: " << sys_realtime() << " sec; CPU: " << sys_cputime() << " sec\n";

	// Return the stderr buffer to its original state.
	std::cerr.rdbuf(old);
	outFile.close();
}
