// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef CONSENSUS_ALIGN_H
#define CONSENSUS_ALIGN_H

#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/consensus.h>

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    char * multipleSequenceAlignment(char * sequences[], char * qualities[], int count,
                                     int bandwidth, int matchScore, int mismatchScore,
                                     int gapOpenScore, int gapExtensionScore);
}


char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities,
                       char oneBaseVsOneGapQualityThreshold);

double getAlignmentIdentity(std::string & seq1, std::string & seq2, int seq1StartPos,
                            int seq1EndPos);

void fillOutQualities(std::vector<std::string> & sequences, std::vector<std::string> & qualities);

void cArrayToCppVector(char * seqArray[], char * qualArray[], int count,
                       std::vector<std::string> & seqVector, std::vector<std::string> & qualVector);

#endif // CONSENSUS_ALIGN_H
