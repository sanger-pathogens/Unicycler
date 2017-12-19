// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "scrub.h"

#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>


char * splitSequences(char * alignmentsString, int seqLength, double startingScore,
                      int posScoreFeatherSize, int negScoreFeatherSize,
                      double posScoreScalingFactor, int splitAdjustment) {

    // Turn the alignments string into a vector of PafAlignment objects.
    std::vector<std::string> alignmentStrings = split(alignmentsString, ';');
    std::vector<PafAlignment> alignments;
    for (auto alignmentString : alignmentStrings)
        alignments.push_back(PafAlignment(alignmentString));

    // Each position in the sequence is scored based on the alignments around it.
    std::vector<double> scores(size_t(seqLength), startingScore);
    std::vector<double> startOverhangScores(size_t(seqLength), 0.0);
    std::vector<double> endOverhangScores(size_t(seqLength), 0.0);

    for (auto a : alignments) {

        // Alignments which span the position earn it positive points (because they suggest
        // contiguity). The larger the span (up to a maximum), the more points are earned.
        int aStart = a.ref_start;
        int aEnd = a.ref_end;
        for (int pos = aStart; pos < aEnd; ++pos) {
            int distanceFromEnd = std::min(pos - aStart, aEnd - pos);
            int relevantDistance = std::min(distanceFromEnd, posScoreFeatherSize);
            double score = posScoreScalingFactor * (double(relevantDistance) / posScoreFeatherSize);
            scores[pos] += score;
        }

        // Tally up the score for start overlaps. Scores are larger for positions in larger
        // overlaps and positions close to where the overlap begins.
        int startOverhangSize = std::min(a.getStartOverhang(), negScoreFeatherSize);
        if (startOverhangSize > 0) {
            double startOverhangRelSize = double(startOverhangSize) / negScoreFeatherSize;
            int rangeStart = aStart - startOverhangSize;
            int rangeEnd = aStart + splitAdjustment;
            for (int pos = rangeStart; pos < rangeEnd; ++pos) {
                int distanceFromClip = aStart - pos;
                double score = double(startOverhangSize - distanceFromClip) / startOverhangSize;
                score = std::min(score, 1.0) * startOverhangRelSize;
                if (pos >= 0 && pos < seqLength)
                    startOverhangScores[pos] -= score;
            }
        }

        // Do the same for end overlaps.
        int endOverhangSize = std::min(a.getEndOverhang(), negScoreFeatherSize);
        if (endOverhangSize > 0) {
            double endOverhangRelSize = double(endOverhangSize) / negScoreFeatherSize;
            int rangeStart = aEnd - splitAdjustment;
            int rangeEnd = aEnd + endOverhangSize;
            for (int pos = rangeStart; pos < rangeEnd; ++pos) {
                int distanceFromClip = pos - aEnd;
                double score = double(endOverhangSize - distanceFromClip) / endOverhangSize;
                score = std::min(score, 1.0) * endOverhangRelSize;
                if (pos >= 0 && pos < seqLength)
                    endOverhangScores[pos] -= score;
            }
        }
    }

    // The final score comes from the positive contribution of spanning alignments and the
    // negative contribution of overhangs. But there needs to be both a start overhang and
    // end overhang to really matter.
    for (int pos = 0; pos < seqLength; ++pos) {
        double overhangScore = std::max(startOverhangScores[pos], endOverhangScores[pos]);
        scores[pos] += overhangScore;
    }

    // Now we get the positive and negative scoring regions of the sequence.
    std::vector<int> positiveScoreRangeStarts, positiveScoreRangeEnds;
    std::vector<int> negativeScoreRangeStarts, negativeScoreRangeEnds;
    int positiveRangeStart = 0;
    bool inPositiveRange = true;
    int negativeRangeStart = 0;
    bool inNegativeRange = true;
    for (int pos = 0; pos < seqLength; ++pos) {
        if (scores[pos] >= 0.0) {
            if (!inPositiveRange) {
                inPositiveRange = true;
                positiveRangeStart = pos;
            }
            if (inNegativeRange) {
                inNegativeRange = false;
                if (pos > negativeRangeStart) {
                    negativeScoreRangeStarts.push_back(negativeRangeStart);
                    negativeScoreRangeEnds.push_back(pos);
                }
            }
        }
        else {  // negative score
            if (inPositiveRange) {
                inPositiveRange = false;
                if (pos > positiveRangeStart) {
                    positiveScoreRangeStarts.push_back(positiveRangeStart);
                    positiveScoreRangeEnds.push_back(pos);
                }
                if (!inNegativeRange) {
                    inNegativeRange = true;
                    negativeRangeStart = pos;
                }
            }
        }
    }
    if (inPositiveRange) {
        positiveScoreRangeStarts.push_back(positiveRangeStart);
        positiveScoreRangeEnds.push_back(seqLength);
    }
    if (inNegativeRange) {
        negativeScoreRangeStarts.push_back(negativeRangeStart);
        negativeScoreRangeEnds.push_back(seqLength);
    }

    // Turn the positive and negative ranges into a string for passing back to Python.
    std::string returnString;
    for (size_t i = 0; i < positiveScoreRangeStarts.size(); ++i) {
        returnString += std::to_string(positiveScoreRangeStarts[i]);
        returnString += "-";
        returnString += std::to_string(positiveScoreRangeEnds[i]);
        if (i < positiveScoreRangeStarts.size() - 1)
            returnString += ",";
    }
    returnString += ";";
    for (size_t i = 0; i < negativeScoreRangeStarts.size(); ++i) {
        returnString += std::to_string(negativeScoreRangeStarts[i]);
        returnString += "-";
        returnString += std::to_string(negativeScoreRangeEnds[i]);
        if (i < negativeScoreRangeStarts.size() - 1)
            returnString += ",";
    }

    return cppStringToCString(returnString);
}


PafAlignment::PafAlignment(std::string alignmentString) {
    std::vector<std::string> parts = split(alignmentString, '\t');
    int read_length = std::stoi(parts[0]);
    read_start = std::stoi(parts[1]);
    read_end = std::stoi(parts[2]);
    int ref_length = std::stoi(parts[3]);
    ref_start = std::stoi(parts[4]);
    ref_end = std::stoi(parts[5]);
    read_end_gap = read_length - read_end;
    ref_end_gap = ref_length - ref_end;
}


int PafAlignment::getStartOverhang() {
    return std::min(read_start, ref_start);
}


int PafAlignment::getEndOverhang() {
    return std::min(read_end_gap, ref_end_gap);
}


// https://stackoverflow.com/questions/236129/split-a-string-in-c
template<typename Out> void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}
