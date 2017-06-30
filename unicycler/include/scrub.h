// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef UNICYCLER_SCRUB_H
#define UNICYCLER_SCRUB_H

#include <string>
#include "string_functions.h"


// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * splitSequences(char * alignmentsString, int seqLength, double startingScore,
                          int posScoreFeatherSize, int negScoreFeatherSize,
                          double posScoreScalingFactor, int splitAdjustment);
}

class PafAlignment {
public:
    PafAlignment(std::string alignmentString);

    int read_start;
    int read_end;
    int ref_start;
    int ref_end;
    int read_end_gap;
    int ref_end_gap;

    int getStartOverhang();
    int getEndOverhang();
};

template<typename Out> void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);

#endif //UNICYCLER_SCRUB_H
