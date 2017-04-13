// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "minimap_align.h"

#include <assert.h>
#include <zlib.h>
#include <iostream>
#include <sstream>

#pragma GCC diagnostic ignored "-Wunused-function"

#include "string_functions.h"
#include "settings.h"

KSEQ_INIT(gzFile, gzread)


char * minimapAlignReads(char * referenceFasta, char * readsFastq, int n_threads,
                         int sensitivityLevel, int preset) {
    // The k-mer size depends on the sensitivity level.
    int k = LEVEL_0_MINIMAP_KMER_SIZE;
    if (sensitivityLevel == 1)
        k = LEVEL_1_MINIMAP_KMER_SIZE;
    else if (sensitivityLevel == 2)
        k = LEVEL_2_MINIMAP_KMER_SIZE;
    else if (sensitivityLevel == 3)
        k = LEVEL_3_MINIMAP_KMER_SIZE;

    // Set up some options and parameters.
    int w = int(.6666667 * k + .499);  // 2/3 of k
    mm_verbose = 0;
    mm_mapopt_t opt;
    mm_mapopt_init(&opt);
	int tbatch_size = 100000000;
	uint64_t ibatch_size = 4000000000ULL;
	float f = 0.001;

    // preset of 0 is default settings.

    // preset of 1 is for mapping reads against themselves: -Sw5 -L100 -m0
    if (preset == 1) {
        opt.flag |= MM_F_AVA | MM_F_NO_SELF;
        opt.min_match = 100;
        opt.merge_frac = 0.0;
        w = 5;
    }
    // preset of 2 is for finding contigs in the string graph: -w5 -L100 -m0
    else if (preset == 2) {
        opt.min_match = 100;
        opt.merge_frac = 0.0;
        w = 5;
    }

    // Redirect minimap's output to a stringstream, instead of outputting it to stdout.
    // http://stackoverflow.com/questions/5419356/redirect-stdout-stderr-to-a-string
    std::stringstream outputBuffer;
    std::streambuf * old = std::cout.rdbuf(outputBuffer.rdbuf());

	bseq_file_t *fp = bseq_open(referenceFasta);
	for (;;) {
		mm_idx_t *mi = 0;
		if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, w, k, MM_IDX_DEF_B, tbatch_size, n_threads, ibatch_size, 1);
		if (mi == 0)
		    break;
		mm_idx_set_max_occ(mi, f);
		mm_map_file(mi, readsFastq, &opt, n_threads, tbatch_size);
		mm_idx_destroy(mi);
	}
	bseq_close(fp);

	// Return the stdout buffer to its original state.
	std::cout.rdbuf(old);

    return cppStringToCString(outputBuffer.str());
}
