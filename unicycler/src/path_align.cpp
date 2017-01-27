// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "path_align.h"

#include <seqan/align.h>
#include "semi_global_align.h"


char * pathAlignment(char * s1, char * s2,
                     int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                     bool useBanding, int bandSize) {

    // Change the sequences to C++ strings.
    std::string sequence1(s1);
    std::string sequence2(s2);

    ScoredAlignment * alignment = pathAlignment(sequence1, sequence2,
                                                matchScore, mismatchScore, gapOpenScore, gapExtensionScore,
                                                useBanding, bandSize);

    if (alignment != 0) {
        std::string returnString = alignment->getFullString();
        delete alignment;
        return cppStringToCString(returnString);
    }
    else
        return cppStringToCString("");
}

// This function runs a mostly-global alignment between two sequences. The only free gaps are those
// at the end of sequence 2.
// It is intended to align a partial path sequence (s1) to a consensus read sequence (s2).
ScoredAlignment * pathAlignment(std::string s1, std::string s2,
                                int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                                bool useBanding, int bandSize) {
    long long startTime = getTime();

    Dna5String sequenceH(s1);
    Dna5String sequenceV(s2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    AlignConfig<false, false, true, false> alignConfig;
    int score;
    if (useBanding) {
        int lowerDiagonal = -bandSize;
        int upperDiagonal = bandSize;
        int lengthDifference = length(s2) - length(s1);

        // If s1 is longer, then we need to expand the upper diagonal a bit.
        if (lengthDifference < 0)
            upperDiagonal -= lengthDifference;

        try {
            score = globalAlignment(alignment, scoringScheme, alignConfig, lowerDiagonal, upperDiagonal);
        }
        catch (...) {
            return 0;
        }
    }
    else {
        try {
            score = globalAlignment(alignment, scoringScheme, alignConfig);
        }
        catch (...) {
            return 0;
        }
    }

    // If the score is too ridiculously low, then something went wrong.
    if (score < -1000000)
        return 0;

    std::string s1Name = "s1";
    std::string s2Name = "s2";

    return new ScoredAlignment(alignment, s1Name, s2Name, s1.length(), s2.length(),
                               0, startTime, 0, true, true, false, scoringScheme);
}


