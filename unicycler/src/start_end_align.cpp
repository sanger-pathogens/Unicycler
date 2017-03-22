// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU Genral Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "start_end_align.h"
#include <seqan/align.h>


// In a start alignment, we expect to find s1 at the start of s2.
int startAlignment(char * s1, char * s2, int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    return startEndAlignment(s1, s2, true, matchScore, mismatchScore, gapOpenScore, gapExtensionScore);
}


// In an end alignment, we expect to find s1 at the end of s2.
int endAlignment(char * s1, char * s2, int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    return startEndAlignment(s1, s2, false, matchScore, mismatchScore, gapOpenScore, gapExtensionScore);
}


int startEndAlignment(char * s1, char * s2, bool start,
                      int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::string sequence1(s1);
    std::string sequence2(s2);

    // Trim the sequence 2 down to a bit more than the expected alignment (to help save time).
    int trim_size = int(sequence1.length() * 1.5);
    int trim_offset;
    if (start)
        trim_offset = 0;
    else
        trim_offset = sequence2.length() - trim_size;
    sequence2 = sequence2.substr(trim_offset, trim_size);

    Dna5String sequenceH(sequence1);
    Dna5String sequenceV(sequence2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    try {
        if (start) {
            AlignConfig<false, false, true, false> alignConfig;  // free end gaps for s1
            globalAlignment(alignment, scoringScheme, alignConfig);
        }
        else {
            AlignConfig<false, true, false, false> alignConfig;  // free end gaps for s2
            globalAlignment(alignment, scoringScheme, alignConfig);
        }
    }
    catch (...) {
        return -1;
    }

    std::ostringstream stream1;
    stream1 << row(alignment, 0);
    std::string seq1Alignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(alignment, 1);
    std::string seq2Alignment =  stream2.str();

    int alignmentLength = std::max(seq1Alignment.size(), seq2Alignment.size());
    if (alignmentLength == 0)
        return -1;

    int seq1Pos = 0, seq2Pos = 0;
    int seq2PosAtSeq1Start = -1, seq2PosAtSeq1End = -1;

    for (int i = 0; i < alignmentLength; ++i) {
        char base1 = seq1Alignment[i];
        char base2 = seq2Alignment[i];

        if (base1 != '-') {
            if (seq2PosAtSeq1Start == -1)
                seq2PosAtSeq1Start = seq2Pos;
            seq2PosAtSeq1End = seq2Pos + 1;
        }

        if (base1 != '-')
            ++seq1Pos;
        if (base2 != '-')
            ++seq2Pos;
    }
    if (start)
        return seq2PosAtSeq1End;
    else
        return seq2PosAtSeq1Start + trim_offset;
}
