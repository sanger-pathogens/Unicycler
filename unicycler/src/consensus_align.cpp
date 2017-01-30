// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "consensus_align.h"

#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include "semi_global_align.h"
#include "string_functions.h"

using namespace seqan;

char * multipleSequenceAlignment(char * sequences[], char * qualities[], int count, int bandwidth,
                                 int matchScore, int mismatchScore, int gapOpenScore,
                                 int gapExtensionScore) {

    // Convert the inputs (arrays of C strings) to C++ vectors, and ensure that the qualities have
    // the same length as their corresponding sequences.
    std::vector<std::string> ungappedSequences, ungappedQualities;
    cArrayToCppVector(sequences, qualities, count, ungappedSequences, ungappedQualities);

    // These vectors will hold the final aligned sequences and qualities.
    std::vector<std::string> gappedSequences, gappedQualities;
    gappedSequences.reserve(count);
    gappedQualities.reserve(count);

    // We want to use a scoring scheme which is almost like the simple scoring scheme we use for
    // read alignment, but with the addition of free Ns.
    typedef Score<int, ScoreMatrix<Dna5, Default> > TScore;
    TScore scoringScheme(gapExtensionScore, gapOpenScore);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            int score;
            if (i == 4 || j == 4) // either base is N
                score = 0;
            else if (i == j)
                score = matchScore;
            else
                score = mismatchScore;
            setScore(scoringScheme, Dna5(i), Dna5(j), score);
        }
    }

    // Prepare data structures for the alignment.
    Align<Dna5String> align;
    resize(rows(align), count);
    for (int i = 0; i < count; ++i)
        assignSource(row(align, i), ungappedSequences[i]);
    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    TStringSet sequenceSet = stringSet(align);
    Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
    String<String<char> > sequenceNames;
    resize(sequenceNames, length(sequenceSet), String<char>("tmpName"));

    MsaOptions<AminoAcid, TScore> msaOpt;
    msaOpt.sc = scoringScheme;
    msaOpt.isDefaultPairwiseAlignment = false;
    msaOpt.method = 0;  // global alignment
    msaOpt.pairwiseAlignmentMethod = 2; // banded
    msaOpt.bandWidth = (unsigned int)bandwidth;

    globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
    convertAlignment(gAlign, align);

    for (int i = 0; i < count; ++i) {
        std::ostringstream stream;
        stream << row(align, i);
        gappedSequences.push_back(stream.str());
    }

    // Add gaps to the quality scores so they match up with the bases.
    int alignmentLength = gappedSequences[0].length();
    for (int i = 0; i < count; ++i) {
        std::string gappedQuality;
        gappedQuality.resize(gappedSequences[i].length(), ' ');
        int pos = 0;
        for (int j = 0; j < alignmentLength; ++j) {
            if (gappedSequences[i][j] != '-')
                gappedQuality[j] = ungappedQualities[i][pos++];
        }
        gappedQualities.push_back(gappedQuality);
    }

    // Before we build a consensus sequence, we may need to deal with one base vs one gap cases,
    // e.g. 'A' vs '-' or '-' vs 'T'. We can't always keep the base, as this would lead to an
    // inflated consensus. So we want a good base quality score threshold, above which we keep the
    // base and below which we keep the gap. This is only necessary if there are places in the
    // alignment with exactly two sequences.
    char oneBaseVsOneGapQualityThreshold = '+';
    if (count == 2) {
        std::vector<char> oneBaseVsOneGapQualities;
        for (int i = 0; i < alignmentLength; ++i) {
            std::vector<char> bases;
            std::vector<char> qualities;
            for (int j = 0; j < count; ++j) {
                char base = toupper(gappedSequences[j][i]);
                char quality = gappedQualities[j][i];
                if (base != 'N') {
                    bases.push_back(base);
                    qualities.push_back(quality);
                }
            }
            if (bases.size() == 2) {
                bool base0IsGap = bases[0] == '-';
                bool base1IsGap = bases[1] == '-';
                if (base0IsGap && !base1IsGap)
                    oneBaseVsOneGapQualities.push_back(qualities[1]);
                else if (base1IsGap && !base0IsGap)
                    oneBaseVsOneGapQualities.push_back(qualities[0]);
            }
        }

        // Set the threshold to the median quality.
        size_t size = oneBaseVsOneGapQualities.size();
        if (size > 0) {
            std::sort(oneBaseVsOneGapQualities.begin(), oneBaseVsOneGapQualities.end());
            if (size % 2 == 0)
                oneBaseVsOneGapQualityThreshold = (oneBaseVsOneGapQualities[size / 2 - 1] +
                                                   oneBaseVsOneGapQualities[size / 2]) / 2;
            else 
                oneBaseVsOneGapQualityThreshold = oneBaseVsOneGapQualities[size / 2];
        }
    }

    // Build a consensus sequence. Sequences are ignored before their first non-N base was seen
    // (for end-only sequences) and after their last non-N base was seen (for start-only
    // sequences).
    std::string consensus, gappedConsensus;
    for (int i = 0; i < alignmentLength; ++i) {
        std::vector<char> bases;
        std::vector<char> qualities;
        bases.reserve(count);
        qualities.reserve(count);
        
        for (int j = 0; j < count; ++j) {
            char base = toupper(gappedSequences[j][i]);
            char quality = gappedQualities[j][i];
            if (base != 'N') {
                bases.push_back(base);
                qualities.push_back(quality);
            }
        }
        if (bases.size() > 0) {
            char mostCommonBase = getMostCommonBase(bases, qualities, oneBaseVsOneGapQualityThreshold);
            if (mostCommonBase != '-')
                consensus.push_back(mostCommonBase);
            gappedConsensus.push_back(mostCommonBase);
        }
        else
            gappedConsensus.push_back('-');
    }

    // Score each sequence against the consensus.
    size_t consensusFirstNonNPos = gappedConsensus.find_first_of("ACGTacgt");
    size_t consensusLastNonNPos = gappedConsensus.find_last_of("ACGTacgt");
    std::vector<double> percentIdentitiesWithConsensus;
    for (int i = 0; i < count; ++i) {
        double identity = getAlignmentIdentity(gappedConsensus, gappedSequences[i],
                                               consensusFirstNonNPos, consensusLastNonNPos);
        percentIdentitiesWithConsensus.push_back(identity);
    }

    std::string returnString = consensus;

    returnString += ';';
    returnString += std::to_string(percentIdentitiesWithConsensus[0]);
    for (int i = 1; i < count; ++i)
        returnString += ',' + std::to_string(percentIdentitiesWithConsensus[i]);

    return cppStringToCString(returnString);
}

char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities, char oneBaseVsOneGapQualityThreshold) {
    std::string baseValues = "ACGT-";

    // Check for the special case of one base vs one gap.
    if (bases.size() == 2) {
        bool base0IsGap = bases[0] == '-';
        bool base1IsGap = bases[1] == '-';

        // If there is one base and one gap and the base has a quality score equal to the
        // threshold, then we choose based on position. I.e. if the quality equals the threshold
        // and the base came first, we return the base, but if it came second, we return the gap.
        if (!base0IsGap && base1IsGap) {
            if (qualities[0] >= oneBaseVsOneGapQualityThreshold) // >= so ties go to the base
                return bases[0]; // base
            else
                return bases[1]; // gap
        }
        else if (base0IsGap && !base1IsGap) {
            if (qualities[1] > oneBaseVsOneGapQualityThreshold) // > so ties go to the gap
                return bases[1]; // base
            else
                return bases[0]; // gap
        }
    }

    // Tally the count for each base.
    std::map<char, int> baseCounts;
    for (int i = 0; i < 5; ++i)
        baseCounts[baseValues[i]] = 0;
    for (size_t i = 0; i < bases.size(); ++i)
        baseCounts[bases[i]]++;

    // If only one base (or a gap) is the most common, return that.
    int largestCount = 0;
    for (int i = 0; i < 5; ++i)
        largestCount = std::max(largestCount, baseCounts[baseValues[i]]);
    std::vector<char> mostCommonBases;
    for (int i = 0; i < 5; ++i) {
        if (baseCounts[baseValues[i]] == largestCount)
            mostCommonBases.push_back(baseValues[i]);
    }
    if (mostCommonBases.size() == 1)
        return mostCommonBases[0];

    // If the code got here, then there's a tie. Check to see if one is a gap, as that should
    // always lose to an actual base.
    if (mostCommonBases.size() == 2 and mostCommonBases[1] == '-')
        return mostCommonBases[0];

    // If the code got here, then there's a tie between two or more bases, so we'll need to use the
    // qualities to decide. The base with the largest Phred sum wins. Since we're just summing the
    // Phred scores, there's no need to worry about +33 vs +64 offset.
    std::map<char, int> phredSums;
    for (int i = 0; i < 4; ++i)
        phredSums[baseValues[i]] = 0;
    for (size_t i = 0; i < bases.size(); ++i)
        phredSums[bases[i]] += int(qualities[i]);
    int largestPhredSum = 0;
    for (size_t i = 0; i < mostCommonBases.size(); ++i)
        largestPhredSum = std::max(largestPhredSum, phredSums[mostCommonBases[i]]);

    // While looping through the bases in their original order, return the first one which has the
    // largest Phred sum and was one of the most common bases. This ensure that a tie (when there
    // are two equally common bases with the same Phred sum) it is broken by the original order in
    // the bases vector, which is arbitrary but at least consistent and not biased towards a
    // particular base.
    for (size_t i = 0; i < bases.size(); ++i) {
        char base = bases[i];
        if (base != '-' && phredSums[base] == largestPhredSum &&
            std::find(mostCommonBases.begin(), mostCommonBases.end(), base) != mostCommonBases.end())
            return base;
    }

    // The code should never get here, as the most common base with the biggest Phred sum should
    // found above.
    return '-';
}


double getAlignmentIdentity(std::string & seq1, std::string & seq2, int seq1StartPos,
                            int seq1EndPos) {
    int alignedLength = 0, matches = 0;

    if (seq1StartPos > seq1EndPos)
        return 0.0;

    for (int i = seq1StartPos; i <= seq1EndPos; ++i) {
        char base1 = seq1[i];
        char base2 = seq2[i];
        if (base1 != '-' || base2 != '-') {
            ++alignedLength;
            if (base1 == base2)
                ++matches;
        }
    }

    return double(matches) / double(alignedLength);
}


// Given a list of sequences and qualities, this function fills out any missing qualities with '+'
// (PHRED+33 for 10% chance of error).
void fillOutQualities(std::vector<std::string> & sequences, std::vector<std::string> & qualities) {
    for (size_t i = 0; i < sequences.size(); ++i)
        qualities[i].resize(sequences[i].length(), '+');
}


void cArrayToCppVector(char * seqArray[], char * qualArray[], int count,
                       std::vector<std::string> & seqVector, std::vector<std::string> & qualVector) {
    seqVector.reserve(count);
    qualVector.reserve(count);
    for (int i = 0; i < count; ++i)
        seqVector.push_back(std::string(seqArray[i]));
    for (int i = 0; i < count; ++i)
        qualVector.push_back(std::string(qualArray[i]));
    fillOutQualities(seqVector, qualVector);
}
