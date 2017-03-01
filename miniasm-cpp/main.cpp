#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include "sys.h"
#include "sdict.h"
#include "miniasm.h"
#include "tclap/CmdLine.h"


#define VERSION "0.3-unicycler"

using namespace TCLAP;
using namespace std;

static void print_subs(const sdict_t *d, const ma_sub_t *sub)
{
	uint32_t i;
	for (i = 0; i < d->n_seq; ++i)
		if (!d->seq[i].del && sub[i].s != sub[i].e)
			printf("%s\t%d\t%d\n", d->seq[i].name, sub[i].s, sub[i].e);
}

static void print_hits(size_t n_hits, const ma_hit_t *hit, const sdict_t *d, const ma_sub_t *sub)
{
	size_t i;
	for (i = 0; i < n_hits; ++i) {
		const ma_hit_t *p = &hit[i];
		const ma_sub_t *rq = &sub[p->qns>>32], *rt = &sub[p->tn];
		printf("%s:%d-%d\t%d\t%d\t%d\t%c\t%s:%d-%d\t%d\t%d\t%d\t%d\t%d\t255\n", d->seq[p->qns>>32].name, rq->s + 1, rq->e, rq->e - rq->s, (uint32_t)p->qns, p->qe,
				"+-"[p->rev], d->seq[p->tn].name, rt->s + 1, rt->e, rt->e - rt->s, p->ts, p->te, p->ml, p->bl);
	}
}

int main(int argc, char *argv[])
{
    string paf_filename, reads_filename, outdir;
    int min_match, min_span, min_dp, min_ovlp, max_hang, gap_fuzz, bub_dist, max_ext, n_rounds;
    float min_iden, int_frac, min_drop, max_drop, final_drop;
    bool prefilter_contained, bi_dir, no_first, no_second;
    try {
        CmdLine cmd("miniasm for Unicycler", ' ', "");

        ValueArg<string> mappingArg(        "p", "paf",        "Input PAF mapping",                                     true,  "",        "paf");
        ValueArg<string> readsArg(          "f", "reads",      "Input long reads",                                      true,  "",        "reads");
        ValueArg<string> outArg(            "",  "out",        "Output directory",                                      false, "output",  "dir");
        SwitchArg prefilterSwitch(          "R", "prefilter",  "prefilter clearly contained reads (2-pass required)",   false);
        ValueArg<int> minMatchLengthArg(    "m", "min_match",  "min match length",                                      false,   100,     "INT");
        ValueArg<float> minIdentityArg(     "i", "min_iden",   "min identity",                                          false,   .05,     "FLOAT");
        ValueArg<int> minSpanArg(           "s", "min_span",   "min span",                                              false,  2000,     "INT");
        ValueArg<int> minCoverageArg(       "c", "min_dp",     "min coverage",                                          false,     1,     "INT");
        ValueArg<int> minOverlapArg(        "o", "min_ovlp",   "min overlap",                                           false,  2000,     "INT");
        ValueArg<int> maxOverHangArg(       "x", "max_hang",   "max over hang length",                                  false,  1000,     "INT");
        ValueArg<float> minEndToEndArg(     "I", "int_frac",   "min end-to-end match ratio",                            false,   0.8,     "FLOAT");
        ValueArg<int> gapFuzzArg(           "g", "gap_fuzz",   "max gap differences between reads for trans-reduction", false,  1000,     "INT");
        ValueArg<int> maxBubbleArg(         "d", "bub_dist",   "max distance for bubble popping",                       false, 50000,     "INT");
        ValueArg<int> smallUnitigArg(       "e", "max_ext",    "small unitig threshold",                                false,     4,     "INT");
        ValueArg<int> overlapRoundsArg(     "n", "n_rounds",   "rounds of short overlap removal",                       false,     2,     "INT");
        ValueArg<float> minOverlapDropArg(  "q", "min_drop",   "min overlap drop ratio",                                false,   0.5,     "FLOAT");
        ValueArg<float> maxOverlapDropArg(  "r", "max_drop",   "max overlap drop ratio",                                false,   0.7,     "FLOAT");
        ValueArg<float> finalOverlapDropArg("F", "final_drop", "aggressive overlap drop ratio in the end",              false,   0.8,     "FLOAT");
        SwitchArg bidirectionalSwitch(      "b", "bi_dir",     "both directions of an arc are present in input",        false);
        SwitchArg skip1PassSwitch(          "1", "no_first",   "skip 1-pass read selection",                            false);
        SwitchArg skip2PassSwitch(          "2", "no_second",  "skip 2-pass read selection",                            false);

        cmd.add(mappingArg);
        cmd.add(readsArg);
        cmd.add(outArg);
        cmd.add(prefilterSwitch);
        cmd.add(minMatchLengthArg);
        cmd.add(minIdentityArg);
        cmd.add(minSpanArg);
        cmd.add(minCoverageArg);
        cmd.add(minOverlapArg);
        cmd.add(maxOverHangArg);
        cmd.add(minEndToEndArg);
        cmd.add(gapFuzzArg);
        cmd.add(maxBubbleArg);
        cmd.add(smallUnitigArg);
        cmd.add(overlapRoundsArg);
        cmd.add(minOverlapDropArg);
        cmd.add(maxOverlapDropArg);
        cmd.add(finalOverlapDropArg);
        cmd.add(bidirectionalSwitch);
        cmd.add(skip1PassSwitch);
        cmd.add(skip2PassSwitch);

        cmd.parse(argc, argv);

        paf_filename = mappingArg.getValue();
        reads_filename = readsArg.getValue();
        outdir = outArg.getValue();
        prefilter_contained = prefilterSwitch.getValue();
        min_match = minMatchLengthArg.getValue();
        min_iden = minIdentityArg.getValue();
        min_span = minSpanArg.getValue();
        min_dp = minCoverageArg.getValue();
        min_ovlp = minOverlapArg.getValue();
        max_hang = maxOverHangArg.getValue();
        int_frac = minEndToEndArg.getValue();
        gap_fuzz = gapFuzzArg.getValue();
        bub_dist = maxBubbleArg.getValue();
        max_ext = smallUnitigArg.getValue();
        n_rounds = overlapRoundsArg.getValue();
        min_drop = minOverlapDropArg.getValue();
        max_drop = maxOverlapDropArg.getValue();
        final_drop = finalOverlapDropArg.getValue();
        bi_dir = !bidirectionalSwitch.getValue();
        no_first = skip1PassSwitch.getValue();
        no_second = skip2PassSwitch.getValue();

    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

    float cov = 40.0;
    string outfmt = "ug";

    mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    string raw_string_graph = outdir + "/01_raw_string_graph.gfa";
    string transitive_reduction_string_graph = outdir + "/02_after_transitive_reduction.gfa";
    string tip_cut_string_graph = outdir + "/03_after_tip_cutting.gfa";
    string bubble_pop_string_graph = outdir + "/04_after_bubble_popping.gfa";
    string cut_overlaps_string_graph_1 = outdir + "/05_after_cutting_short_overlaps.gfa";
    string remove_internal_string_graph = outdir + "/06_after_removing_internal_seqs.gfa";
    string cut_overlaps_string_graph_2 = outdir + "/07_after_cutting_short_overlaps.gfa";
    string final_string_graph = outdir + "/08_final_string_graph.gfa";
    string final_unitig_graph = outdir + "/09_unitig_graph.gfa";

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
    ma_hit_t * hits = read_hits_file(paf_filename.c_str(), min_span, min_match, read_dict, &num_hits, bi_dir, excluded_reads);
    std::cerr << "\n";

    ma_sub_t *subreads = 0;
	if (!no_first) {
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
	}

	if (!no_second) {
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

        // Toss out chimeric reads.
		remove_chimeric_reads(max_hang, min_dp, num_hits, hits, read_dict, subreads);

        // Toss out contained reads (this is a big one and gets rid of a lot).
		num_hits = remove_contained_reads(max_hang, int_frac, min_ovlp, read_dict, subreads, num_hits, hits);
        std::cerr << "\n";
	}

	hits = (ma_hit_t*)realloc(hits, num_hits * sizeof(ma_hit_t));

//	if (outfmt == "bed") {
//		print_subs(read_dict, subreads);
//	}
//    else if (outfmt == "paf") {
//		print_hits(num_hits, hits, read_dict, subreads);
//	}
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

    ma_ug_t * unitig_graph = 0;
    cerr << "===> Step 5: generating unitigs <===\n";
    unitig_graph = make_unitig_graph(string_graph);
    if (reads_filename.c_str())
        generate_unitig_seqs(unitig_graph, read_dict, subreads, reads_filename.c_str());
    save_unitig_graph(unitig_graph, read_dict, subreads, final_unitig_graph);
    destroy_unitig_graph(unitig_graph);
    std::cerr << "\n";

    save_string_graph(string_graph, read_dict, subreads, final_string_graph, reads_filename.c_str());
    destroy_string_graph(string_graph);

    // Clean up!
	free(subreads); free(hits);
	destroy_seq_dict(read_dict);
	if (excluded_reads)
        destroy_seq_dict(excluded_reads);

    cerr << "Version: " << VERSION << "\n";
    cerr << "CMD:"; for (int i = 0; i < argc; ++i) cerr << " " << argv[i]; cerr << "\n";
    cerr << "Real time: " << sys_realtime() << " sec; CPU: " << sys_cputime() << " sec\n";
	return 0;
}
