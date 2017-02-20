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

#include "string_functions.h"
#include "settings.h"

KSEQ_INIT(gzFile, gzread)


char * minimapAlignReads(char * referenceFasta, char * readsFastq, int n_threads,
                         int sensitivityLevel) {

    int k = LEVEL_0_MINIMAP_KMER_SIZE;
    if (sensitivityLevel == 1)
        k = LEVEL_1_MINIMAP_KMER_SIZE;
    else if (sensitivityLevel == 2)
        k = LEVEL_2_MINIMAP_KMER_SIZE;
    else if (sensitivityLevel == 3)
        k = LEVEL_3_MINIMAP_KMER_SIZE;

    // Load in the reads.
    gzFile f = gzopen(readsFastq, "r");
    assert(f);
    kseq_t *ks = kseq_init(f);

    // Build the reference index.
    int w = (int)(.6666667 * k + .499);  // 2/3 of k
    mm_verbose = 0;
    mm_idx_t * mi = mm_idx_build(referenceFasta, w, k, n_threads);
    assert(mi);

    std::string outputString;

    mm_mapopt_t opt;
    mm_mapopt_init(&opt); // initialize mapping parameters
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        const mm_reg1_t *reg;
        int j, n_reg;
        // get all hits for the query
        reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &opt, 0);
        // traverse hits and print them out
        for (j = 0; j < n_reg; ++j) {
            const mm_reg1_t *r = &reg[j];

            std::string readName = ks->name.s;
            // int readLen = ks->seq.l;
            int readStart = r->qs;
            int readEnd = r->qe;
            char readStrand = "+-"[r->rev];

            std::string refName = mi->name[r->rid];
            // int refLen = mi->len[r->rid];
            int refStart = r->rs;
            int refEnd = r->re;

            int numberOfMinimisers = r->cnt;

            std::string alignmentString = readName + "\t";
            alignmentString += std::to_string(readStart) + "\t" + std::to_string(readEnd) + "\t" + readStrand + "\t";
            alignmentString += refName + "\t";
            alignmentString += std::to_string(refStart) + "\t" + std::to_string(refEnd) + "\t";
            alignmentString += std::to_string(numberOfMinimisers) + "\t";
            alignmentString += "\n";

            outputString += alignmentString;
        }
    }
    mm_tbuf_destroy(tbuf);

    mm_idx_destroy(mi);
    kseq_destroy(ks);
    gzclose(f);

    return cppStringToCString(outputString);
}
