// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "random_alignments.h"

#include <seqan/align.h>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>
#include <random>
#include <thread>

#include "semi_global_align.h"
#include "global_align.h"



// This function runs a bunch of alignments between random sequences to get a mean and std dev of
// the scaled scores. It return them in a C string (for Python).
char * getRandomSequenceAlignmentScores(int seqLength, int n,
                                        int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::vector<double> scores;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 3);

    for (int i = 0; i < n; ++i) {
        std::string s1 = getRandomSequence(seqLength, gen, dist);
        std::string s2 = getRandomSequence(seqLength, gen, dist);
        ScoredAlignment * alignment = fullyGlobalAlignment(s1, s2, matchScore, mismatchScore, gapOpenScore, gapExtensionScore);

        if (alignment != 0) {
            scores.push_back(alignment->m_scaledScore);
            delete alignment;
        }
    }

    double mean = 0.0, stdev = 0.0;
    getMeanAndStDev(scores, mean, stdev);
    return cppStringToCString(std::to_string(mean) + "," + std::to_string(stdev));
}

// This function returns lots of information about random global alignments.
char * getRandomSequenceAlignmentErrorRates(int seqLength, int n,
                                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::vector<double> matchesOverAlignmentLength;
    std::vector<double> errorsOverAlignmentLength;
    std::vector<double> matchesOverSeqLength;
    std::vector<double> errorsOverSeqLength;
    std::vector<double> errorsOverSeqLengthCountAllInsertionsAsOne;
    std::vector<double> scores;
    std::vector<double> alignmentLengths;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 3);

    for (int i = 0; i < n; ++i) {
        int totalMatches = 0;
        int totalErrors = 0;
        int totalErrorsCountAllInsertionsAsOne = 0;

        std::string s1 = getRandomSequence(seqLength, gen, dist);
        std::string s2 = getRandomSequence(seqLength, gen, dist);

        Dna5String sequenceH(s1);
        Dna5String sequenceV(s2);

        Align<Dna5String, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

        AlignConfig<false, false, false, false> alignConfig;
        int score = globalAlignment(alignment, scoringScheme, alignConfig);

        // Extract the alignment sequences into C++ strings for constant time random access.
        std::ostringstream stream1;
        stream1 << row(alignment, 0);
        std::string s1Alignment =  stream1.str();
        std::ostringstream stream2;
        stream2 << row(alignment, 1);
        std::string s2Alignment =  stream2.str();
        int alignmentLength = s1Alignment.size();

        int insertionLength = 0;
        for (size_t i = 0; i < s1Alignment.size(); ++i)
        {
            if (s1Alignment[i] != '-' && s2Alignment[i] == '-')
                ++insertionLength;
            else
                insertionLength = 0;

            if (s1Alignment[i] == '-' || s2Alignment[i] == '-' || s1Alignment[i] != s2Alignment[i]) {
                ++totalErrors;
                if (insertionLength < 2)
                    ++totalErrorsCountAllInsertionsAsOne;
            }
            else // match
                ++totalMatches;
        }
        matchesOverAlignmentLength.push_back(double(totalMatches) / double(alignmentLength));
        errorsOverAlignmentLength.push_back(double(totalErrors) / double(alignmentLength));
        matchesOverSeqLength.push_back(double(totalMatches) / double(seqLength));
        errorsOverSeqLength.push_back(double(totalErrors) / double(seqLength));
        errorsOverSeqLengthCountAllInsertionsAsOne.push_back(double(totalErrorsCountAllInsertionsAsOne) / double(seqLength));
        scores.push_back(double(score));
        alignmentLengths.push_back(double(s1Alignment.size()));
    }

    std::string returnString;
    double mean, stdev;

    returnString += "identity (using alignment length) mean\t";
    returnString += "identity (using alignment length) std dev\t";
    returnString += "errors (using alignment length) mean\t";
    returnString += "errors (using alignment length) std dev\t";
    returnString += "identity (using sequence length) mean\t";
    returnString += "identity (using sequence length) std dev\t";
    returnString += "errors (using sequence length) mean\t";
    returnString += "errors (using sequence length) std dev\t";
    returnString += "errors (using sequence length, count all insertions as one) mean\t";
    returnString += "errors (using sequence length, count all insertions as one) std dev\t";
    returnString += "score mean\t";
    returnString += "score std dev\t";
    returnString += "alignment length mean\t";
    returnString += "alignment length std dev\n";

    getMeanAndStDev(matchesOverAlignmentLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverAlignmentLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(matchesOverSeqLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverSeqLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverSeqLengthCountAllInsertionsAsOne, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(scores, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(alignmentLengths, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\n";

    return cppStringToCString(returnString);
}

std::string getRandomSequence(int seqLength, std::mt19937 & gen, std::uniform_int_distribution<int> & dist) {
    std::string seq;
    seq.reserve(seqLength);
    for (int i = 0 ; i < seqLength; ++i)
        seq += getRandomBase(gen, dist);
    return seq;
}

char getRandomBase(std::mt19937 & gen, std::uniform_int_distribution<int> & dist) {
    int baseNum = dist(gen);
    if (baseNum == 0)
        return 'A';
    else if (baseNum == 1)
        return 'C';
    else if (baseNum == 2)
        return 'G';
    else // baseNum == 3
        return 'T';
}

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdDev) {
    mean = 0.0;
    stdDev = 0.0;
    int count = v.size();
    if (count < 1)
        return;
    double devSum = 0.0;
    for (int i = 0; i < count; ++i)
        mean += v[i];
    mean /= count;
    for (int i = 0; i < count; ++i) {
        double dev = v[i] - mean;
        devSum += dev * dev;
    }
    stdDev = sqrt(devSum / v.size());
}

char * simulateDepths(int alignmentLengths[], int alignmentCount, int refLength, int iterations, int threadCount) {


    std::vector<int> minDepthCounts;
    std::vector<int> maxDepthCounts;
    std::mutex mut;

    std::vector<std::thread *> threads;
    int iterationsPerThread = iterations / threadCount;
    int iterationsInFirstThread = iterations - (iterationsPerThread * (threadCount - 1));
    for (int i = 0; i < threadCount; ++i) {
        int iterationsThisThread;
        if (i == 0)
            iterationsThisThread = iterationsInFirstThread;
        else
            iterationsThisThread = iterationsPerThread;
        std::thread * thread = new std::thread(simulateDepthsOneThread, alignmentLengths, alignmentCount, refLength, iterationsThisThread,
                                               &minDepthCounts, &maxDepthCounts, &mut);
        threads.push_back(thread);
    }
    for (int i = 0; i < threadCount; ++i) {
        threads[i]->join();
        delete threads[i];
    }


    std::vector<double> minDepthDistribution;
    for (size_t i = 0; i < minDepthCounts.size(); ++i)
        minDepthDistribution.push_back(double(minDepthCounts[i]) / iterations);

    std::vector<double> maxDepthDistribution;
    for (size_t i = 0; i < maxDepthCounts.size(); ++i)
        maxDepthDistribution.push_back(double(maxDepthCounts[i]) / iterations);

    std::string returnString;
    for (size_t i = 0; i < minDepthDistribution.size(); ++i) {
        if (i > 0)
            returnString += ',';
        returnString += std::to_string(i) + ':' + toStringWithPrecision(minDepthDistribution[i], 12);
    }
    returnString += ';';
    for (size_t i = 0; i < maxDepthDistribution.size(); ++i) {
        if (i > 0)
            returnString += ',';
        returnString += std::to_string(i) + ':' + toStringWithPrecision(maxDepthDistribution[i], 12);
    }
    return cppStringToCString(returnString);
}

void simulateDepthsOneThread(int alignmentLengths[], int alignmentCount, int refLength, int iterations,
                             std::vector<int> * minDepthCounts, std::vector<int> * maxDepthCounts,
                             std::mutex * mut) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, refLength-1);

    for (int i = 0; i < iterations; ++i) {

        std::vector<std::pair<int, int> > alignmentPositions;

        for (int j = 0; j < alignmentCount; ++j) {
            int alignmentLength = alignmentLengths[j];
            int start = dist(gen);
            int end = start + alignmentLength;

            // If the end landed within the reference, then we just need to make a single alignment.
            if (end <= refLength) {
                alignmentPositions.push_back(std::pair<int, int>(start, 1));
                alignmentPositions.push_back(std::pair<int, int>(end, -1));
            }
            // If the end landed past the reference, then we make two alignments (looping it back to the start).
            else {
                alignmentPositions.push_back(std::pair<int, int>(start, 1));
                alignmentPositions.push_back(std::pair<int, int>(refLength, -1));
                alignmentPositions.push_back(std::pair<int, int>(0, 1));
                alignmentPositions.push_back(std::pair<int, int>(end - refLength, -1));
            }
        }
        std::sort(alignmentPositions.begin(), alignmentPositions.end());

        int lastPos = 0;
        int currentDepth = 0;
        int minDepth = std::numeric_limits<int>::max();
        int maxDepth = std::numeric_limits<int>::min();
        for (size_t j = 0; j < alignmentPositions.size(); ++j)
        {
            int pos = alignmentPositions[j].first;
            int change = alignmentPositions[j].second;

            int basesAtCurrentDepth = pos - lastPos;
            if (basesAtCurrentDepth > 0) {
                minDepth = std::min(minDepth, currentDepth);
                maxDepth = std::max(maxDepth, currentDepth);
            }
            currentDepth += change;
            lastPos = pos;
        }
        if (lastPos != refLength) {
            int basesAtCurrentDepth = refLength - lastPos;
            if (basesAtCurrentDepth > 0) {
                minDepth = std::min(minDepth, currentDepth);
                maxDepth = std::max(maxDepth, currentDepth);
            }
        }

        mut->lock();
        if (int(minDepthCounts->size()) < minDepth + 1)
            minDepthCounts->resize(minDepth + 1, 0);
        (*minDepthCounts)[minDepth] += 1;
        if (int(maxDepthCounts->size()) < maxDepth + 1)
            maxDepthCounts->resize(maxDepth + 1, 0);
        (*maxDepthCounts)[maxDepth] += 1;
        mut->unlock();
    }
}

// http://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
std::string toStringWithPrecision(double val, int precision)
{
    std::ostringstream out;
    out << std::setprecision(precision) << val;
    return out.str();
}











