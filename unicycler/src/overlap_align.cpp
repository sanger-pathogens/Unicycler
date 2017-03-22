// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU Genral Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "overlap_align.h"
#include "string_functions.h"
#include <seqan/align.h>


char * overlapAlignment(char * s1, char * s2,
                        int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                        int guessOverlap) {
    std::string sequence1(s1);
    std::string sequence2(s2);

    // Trim the sequences down to a bit more than the expected overlap. This will help save time
    // because the non-overlapping sequence isn't informative.
    int trim_size = int((guessOverlap + 100) * 1.5);
    if (trim_size < int(sequence1.length()))
        sequence1 = sequence1.substr(sequence1.length() - trim_size, trim_size);
    if (trim_size < int(sequence2.length()))
        sequence2 = sequence2.substr(0, trim_size);

    Dna5String sequenceH(sequence1);
    Dna5String sequenceV(sequence2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The alignment configuration allows for overlap from s1 -> s2.
    AlignConfig<true, false, true, false> alignConfig;
    try {
        globalAlignment(alignment, scoringScheme, alignConfig);
    }
    catch (...) {
        return cppStringToCString("-1,-1");
    }

    std::ostringstream stream1;
    stream1 << row(alignment, 0);
    std::string seq1Alignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(alignment, 1);
    std::string seq2Alignment =  stream2.str();

    int alignmentLength = std::max(seq1Alignment.size(), seq2Alignment.size());
    if (alignmentLength == 0)
        return cppStringToCString("-1,-1");

    int seq1Pos = 0, seq2Pos = 0;
    int seq1PosAtSeq2Start = -1, seq2PosAtSeq1End = -1;

    for (int i = 0; i < alignmentLength; ++i) {
        char base1 = seq1Alignment[i];
        char base2 = seq2Alignment[i];

        if (base1 != '-')
            seq2PosAtSeq1End = seq2Pos + 1;
        if (base2 != '-' && seq1PosAtSeq2Start == -1)
            seq1PosAtSeq2Start = seq1Pos;

        if (base1 != '-')
            ++seq1Pos;
        if (base2 != '-')
            ++seq2Pos;
    }

    int overlap1 = seq1Pos - seq1PosAtSeq2Start;
    int overlap2 = seq2PosAtSeq1End;
    return cppStringToCString(std::to_string(overlap1) + "," + std::to_string(overlap2));
}
