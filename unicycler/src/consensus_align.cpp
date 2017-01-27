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
#include "semi_global_align.h"
#include <algorithm>
#include <cmath>

char * multipleSequenceAlignment(char * fullSpanSequences[], char * fullSpanQualities[], int fullSpanCount, 
                                 char * startOnlySequences[], char * startOnlyQualities[], int startOnlyCount, 
                                 char * endOnlySequences[], char * endOnlyQualities[], int endOnlyCount, 
                                 int bandwidth, int matchScore, int mismatchScore,
                                 int gapOpenScore, int gapExtensionScore) {

    // Convert the inputs (arrays of C strings) to C++ vectors, and ensure that the qualities have
    // the same length as their corresponding sequences.
    std::vector<std::string> fullSpanSeqs, fullSpanQuals, startOnlySeqs, startOnlyQuals, endOnlySeqs, endOnlyQuals;
    cArrayToCppVector(fullSpanSequences, fullSpanQualities, fullSpanCount, fullSpanSeqs, fullSpanQuals);
    cArrayToCppVector(startOnlySequences, startOnlyQualities, startOnlyCount, startOnlySeqs, startOnlyQuals);
    cArrayToCppVector(endOnlySequences, endOnlyQualities, endOnlyCount, endOnlySeqs, endOnlyQuals);

    // Pad out the start/end sequences/qualities with N/+
    int totalFullSpanLength = 0;
    for (int i = 0; i < fullSpanCount; ++i)
        totalFullSpanLength += fullSpanSeqs[i].length();
    int meanFullSpanLength = int(0.5 + (double(totalFullSpanLength) / fullSpanCount));
    padToLength(startOnlySeqs, startOnlyQuals, meanFullSpanLength, false);
    padToLength(endOnlySeqs, endOnlyQuals, meanFullSpanLength, true);

    // Make vectors of all sequences/qualities together.
    int totalSeqCount = fullSpanCount + startOnlyCount + endOnlyCount;
    std::vector<std::string> ungappedSequences, ungappedQualities;
    ungappedSequences.reserve(totalSeqCount);
    ungappedQualities.reserve(totalSeqCount);
    ungappedSequences.insert(ungappedSequences.end(), fullSpanSeqs.begin(), fullSpanSeqs.end());
    ungappedQualities.insert(ungappedQualities.end(), fullSpanQuals.begin(), fullSpanQuals.end());
    ungappedSequences.insert(ungappedSequences.end(), startOnlySeqs.begin(), startOnlySeqs.end());
    ungappedQualities.insert(ungappedQualities.end(), startOnlyQuals.begin(), startOnlyQuals.end());
    ungappedSequences.insert(ungappedSequences.end(), endOnlySeqs.begin(), endOnlySeqs.end());
    ungappedQualities.insert(ungappedQualities.end(), endOnlyQuals.begin(), endOnlyQuals.end());

    // These vectors will hold the final aligned sequences and qualities.
    std::vector<std::string> gappedSequences, gappedQualities;
    gappedSequences.reserve(totalSeqCount);
    gappedQualities.reserve(totalSeqCount);

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
    resize(rows(align), totalSeqCount);
    for (int i = 0; i < totalSeqCount; ++i)
        assignSource(row(align, i), ungappedSequences[i]);
    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    TStringSet sequenceSet = stringSet(align);
    Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
    String<String<char> > sequenceNames;
    resize(sequenceNames, length(sequenceSet), String<char>("tmpName"));
    MsaOptions<AminoAcid, TScore> msaOpt;
    msaOpt.sc = scoringScheme;
    msaOpt.isDefaultPairwiseAlignment = false;
    msaOpt.pairwiseAlignmentMethod = 2;
    msaOpt.bandWidth = bandwidth;

    // This calls my custom copy of the the globalMsaAlignment function. It does global alignments
    // between pairs of full-span sequences (specified by fullLengthCount) and overlap alignments
    // between all other pairs (full-span to partial pairs and partial to partial pairs).
    globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt, fullSpanCount, startOnlyCount, endOnlyCount);

    convertAlignment(gAlign, align);

    // std::cout << "\n" << align << "\n"; // TEMP

    for (int i = 0; i < totalSeqCount; ++i) {
        std::ostringstream stream;
        stream << row(align, i);
        gappedSequences.push_back(stream.str());
    }

    // Add gaps to the quality scores so they match up with the bases.
    int alignmentLength = gappedSequences[0].length();
    for (int i = 0; i < totalSeqCount; ++i) {
        // std::cout << gappedSequences[i] << "\n"; // TEMP
        std::string gappedQuality;
        gappedQuality.resize(gappedSequences[i].length(), ' ');
        int pos = 0;
        for (int j = 0; j < alignmentLength; ++j) {
            if (gappedSequences[i][j] != '-')
                gappedQuality[j] = ungappedQualities[i][pos++];
        }
        gappedQualities.push_back(gappedQuality);
        // std::cout << gappedQualities[i] << "\n"; // TEMP
    }

    // For each gapped sequence, get the position of the first and last non-N base. This is useful
    // for sequences which don't span the whole alignment but are just at the start or end.
    std::vector<int> firstNonN, lastNonN;
    firstNonN.reserve(totalSeqCount);
    lastNonN.reserve(totalSeqCount);
    for (int i = 0; i < totalSeqCount; ++i) {
        size_t firstNonNPos = gappedSequences[i].find_first_of("ACGTacgt");
        size_t lastNonNPos = gappedSequences[i].find_last_of("ACGTacgt");
        if (firstNonNPos != std::string::npos) {
            firstNonN.push_back(int(firstNonNPos));
            lastNonN.push_back(int(lastNonNPos));
        }
        else {
            firstNonN.push_back(0);
            lastNonN.push_back(gappedSequences[i].length() - 1);
        }
        // std::cout << firstNonN.back() << "  " << lastNonN.back() << "\n"; // TEMP
    }

    // The final consensus alignment should begin at the point where the majority of full-span
    // sequences start/end.
    std::vector<int> fullSpanStarts, fullSpanEnds;
    for (int i = 0; i < fullSpanCount; ++i) {
        fullSpanStarts.push_back(firstNonN[i]);
        fullSpanEnds.push_back(lastNonN[i]);
    }
    std::sort(fullSpanStarts.begin(), fullSpanStarts.end());
    int consensusStart = fullSpanStarts[fullSpanCount / 2];
    std::sort(fullSpanEnds.begin(), fullSpanEnds.end(), std::greater<int>());
    int consensusEnd = fullSpanEnds[fullSpanCount / 2];


    // Before we build a consensus sequence, we may need to deal with one base vs one gap cases,
    // e.g. 'A' vs '-' or '-' vs 'T'. We can't always keep the base, as this would lead to an
    // inflated consensus. So we want a good base quality score threshold, above which we keep the
    // base and below which we keep the gap. This is only necessary if there are places in the
    // alignment with exactly two sequences.
    char oneBaseVsOneGapQualityThreshold = '+';
    if ((fullSpanCount == 2) ||
        (fullSpanCount == 1 && (startOnlyCount > 0 || endOnlyCount > 0))) {
        std::vector<char> oneBaseVsOneGapQualities;
        for (int i = 0; i < alignmentLength; ++i) {
            std::vector<char> bases;
            std::vector<char> qualities;
            for (int j = 0; j < totalSeqCount; ++j) {
                if (i < firstNonN[j] || i > lastNonN[j])
                    continue;
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
        // std::cout << "oneBaseVsOneGapQualityThreshold: " << oneBaseVsOneGapQualityThreshold << "\n"; // TEMP
    }

    // Build a consensus sequence. Sequences are ignored before their first non-N base was seen
    // (for end-only sequences) and after their last non-N base was seen (for start-only
    // sequences).
    std::string consensus, gappedConsensus;
    for (int i = 0; i < alignmentLength; ++i) {
        std::vector<char> bases;
        std::vector<char> qualities;
        bases.reserve(totalSeqCount);
        qualities.reserve(totalSeqCount);

        // std::cout << "Position:  " << i << "\n"; // TEMP

        for (int j = 0; j < totalSeqCount; ++j) {
            if (i < firstNonN[j] || i > lastNonN[j])
                continue;
            char base = toupper(gappedSequences[j][i]);
            char quality = gappedQualities[j][i];
            if (base != 'N') {
                bases.push_back(base);
                qualities.push_back(quality);
            }
        }
        if (bases.size() > 0 && i >= consensusStart && i <= consensusEnd) {
            char mostCommonBase = getMostCommonBase(bases, qualities, oneBaseVsOneGapQualityThreshold);
            // std::cout << "Call:      " << mostCommonBase << "\n\n"; // TEMP
            if (mostCommonBase != '-')
                consensus.push_back(mostCommonBase);
            gappedConsensus.push_back(mostCommonBase);
        }
        else
            gappedConsensus.push_back('-');
    }
    // std::cout << "\n" << gappedConsensus << "\n"; // TEMP

    // Score each sequence against the consensus.
    size_t consensusFirstNonNPos = gappedConsensus.find_first_of("ACGTacgt");
    size_t consensusLastNonNPos = gappedConsensus.find_last_of("ACGTacgt");
    std::vector<double> percentIdentitiesWithConsensus;
    for (int i = 0; i < totalSeqCount; ++i) {
        double identity = getAlignmentIdentity(gappedConsensus, gappedSequences[i],
                                               consensusFirstNonNPos, consensusLastNonNPos,
                                               firstNonN[i], lastNonN[i]);
        percentIdentitiesWithConsensus.push_back(identity);
    }

    // // Score each sequence against each other.
    // std::vector< std::vector<double> > percentIdentitiesBetweenReads;
    // for (int i = 0; i < totalSeqCount; ++i)
    //     percentIdentitiesBetweenReads.push_back(std::vector<double>(totalSeqCount, 0.0));
    // for (int i = 0; i < totalSeqCount; ++i) {
    //     for (int j = i; j < totalSeqCount; ++j) {
    //         if (i == j)
    //             percentIdentitiesBetweenReads[i][j] = 1.0;
    //         else {
    //             double identity = getAlignmentIdentity(gappedSequences[i], gappedSequences[j],
    //                                                    firstNonN[i], lastNonN[i], firstNonN[j], lastNonN[j]);
    //             percentIdentitiesBetweenReads[i][j] = identity;
    //             percentIdentitiesBetweenReads[j][i] = identity;
    //         }
    //     }
    // }

    std::string returnString = consensus;

    returnString += ';';
    returnString += std::to_string(percentIdentitiesWithConsensus[0]);
    for (int i = 1; i < totalSeqCount; ++i)
        returnString += ',' + std::to_string(percentIdentitiesWithConsensus[i]);

    // returnString += ';';
    // for (int i = 0; i < totalSeqCount; ++i) {
    //     for (int j = 0; j < totalSeqCount; ++j) {
    //         returnString += std::to_string(percentIdentitiesBetweenReads[i][j]);
    //         if (j != totalSeqCount - 1)
    //             returnString += ',';
    //     }
    //     if (i != totalSeqCount - 1)
    //         returnString += ';';
    // }
    return cppStringToCString(returnString);
}

char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities, char oneBaseVsOneGapQualityThreshold) {
    std::string baseValues = "ACGT-";

    // std::cout << "Bases:     "; // TEMP
    // for (size_t i = 0; i < bases.size(); ++i) // TEMP
    //     std::cout << bases[i]; // TEMP
    // std::cout << "\n"; // TEMP

    // std::cout << "Qualities: "; // TEMP
    // for (size_t i = 0; i < qualities.size(); ++i) // TEMP
    //     std::cout << qualities[i]; // TEMP
    // std::cout << "\n"; // TEMP

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

    // std::cout << "Counts:    "; // TEMP
    // std::cout << "A: " << baseCounts['A'] << "  "; // TEMP
    // std::cout << "C: " << baseCounts['C'] << "  "; // TEMP
    // std::cout << "G: " << baseCounts['G'] << "  "; // TEMP
    // std::cout << "T: " << baseCounts['T'] << "  "; // TEMP
    // std::cout << "-: " << baseCounts['-'] << "\n"; // TEMP

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


double getAlignmentIdentity(std::string & seq1, std::string & seq2,
                            int seq1StartPos, int seq1EndPos, int seq2StartPos, int seq2EndPos) {
    // std::cout << "\n"; // TEMP
    // std::cout << seq1 << "\n"; // TEMP
    // std::cout << seq1StartPos << ", " << seq1EndPos << "\n"; // TEMP
    // std::cout << seq2 << "\n"; // TEMP
    // std::cout << seq2StartPos << ", " << seq2EndPos << "\n"; // TEMP
    // std::cout << "\n"; // TEMP

    int alignedLength = 0, matches = 0;
    int start = std::max(seq1StartPos, seq2StartPos);
    int end = std::min(seq1EndPos, seq2EndPos);

    if (start > end)
        return 0.0;

    for (int i = start; i <= end; ++i) {
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


// If the input string is less than the given length, this function will add Ns to fill it out. If
// it is longer, this function will trim it down.
void padToLength(std::vector<std::string> & sequences, std::vector<std::string> & qualities, int length, bool putAtStart) {
    for (size_t i = 0; i < sequences.size(); ++i) {
        int missingBases = length - sequences[i].length();
        if (missingBases > 0) {
            std::string newBases(missingBases, 'N');
            std::string newQualities(missingBases, '+');
            if (putAtStart) {
                sequences[i] = newBases + sequences[i];
                qualities[i] = newQualities + qualities[i];
            }
            else {
                sequences[i] = sequences[i] + newBases;
                qualities[i] = qualities[i] + newQualities;
            }
        }
        else if (missingBases < 0) {
            if (putAtStart) {
                int charsToRemove = -missingBases;
                sequences[i] = sequences[i].substr(charsToRemove);
                qualities[i] = qualities[i].substr(charsToRemove);
            }
            else {
                sequences[i] = sequences[i].substr(0, length);
                qualities[i] = qualities[i].substr(0, length);
            }
        }
    }
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