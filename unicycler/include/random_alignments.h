// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef RANDOM_ALIGNMENTS_H
#define RANDOM_ALIGNMENTS_H

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "scoredalignment.h"
#include <mutex>

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * getRandomSequenceAlignmentScores(int seqLength, int n,
                                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
    char * getRandomSequenceAlignmentErrorRates(int seqLength, int n,
                                               int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
    char * simulateDepths(int alignmentLengths[], int alignmentCount, int refLength, int iterations, int threadCount);
}

void simulateDepthsOneThread(int alignmentLengths[], int alignmentCount, int refLength, int iterations,
                             std::vector<int> * minDepthCounts, std::vector<int> * maxDepthCounts,
                             std::mutex * mut);

std::string getRandomSequence(int seqLength, std::mt19937 & gen, std::uniform_int_distribution<int> & dist);

char getRandomBase(std::mt19937 & gen, std::uniform_int_distribution<int> & dist);

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);

std::string toStringWithPrecision(double val, int precision);


#endif // RANDOM_ALIGNMENTS_H
