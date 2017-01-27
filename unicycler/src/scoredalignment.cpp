// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "scoredalignment.h"

#include <iostream>

ScoredAlignment::ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment, 
                                 std::string & readName, std::string & refName,
                                 int readLength, int refLength,
                                 int refOffset, long long startTime, int bandSize,
                                 bool startImmediately, bool goToEndSeq1, bool goToEndSeq2,
                                 Score<int, Simple> & scoringScheme):
    m_readName(readName), m_refName(refName), m_readLength(readLength), m_refLength(refLength),
    m_readStartPos(-1), m_refStartPos(-1), m_rawScore(0), m_bandSize(bandSize)
{
    // Extract the alignment sequences into C++ strings for constant time random access.
    std::ostringstream stream1;
    stream1 << row(alignment, 0);
    std::string readAlignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(alignment, 1);
    std::string refAlignment =  stream2.str();

    int alignmentLength = std::max(readAlignment.size(), refAlignment.size());
    if (alignmentLength == 0)
        return;

    CigarType currentCigarType;
    int currentCigarLength = 0;
    int readBases = 0, refBases = 0;
    std::vector<CigarType> cigarTypes;
    std::vector<int> cigarLengths;

    bool alignmentStarted = false;
    bool readStarted = false, refStarted = false;
    int alignmentStartPos = -1;
    int alignmentEndPos = -1;

    if (startImmediately) {
        alignmentStarted = true;
        readStarted = true;
        refStarted = true;
        m_readStartPos = 0;
        m_refStartPos = 0;
        alignmentStartPos = 0;
    }

    for (int i = 0; i < alignmentLength; ++i) {
        char base1 = readAlignment[i];
        char base2 = refAlignment[i];

        // We consider the alignment to have started when we've encountered a base in both
        // sequences (though not necessarily at the same time).
        if (base1 != '-')
            readStarted = true;
        if (base2 != '-')
            refStarted = true;
        if (readStarted && refStarted && !alignmentStarted) {
            m_readStartPos = readBases;
            m_refStartPos = refBases;
            alignmentStarted = true;
            alignmentStartPos = i;
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (i == 0)
            currentCigarType = cigarType;
        if (cigarType == currentCigarType)
            ++currentCigarLength;
        else {
            cigarTypes.push_back(currentCigarType);
            cigarLengths.push_back(currentCigarLength);
            currentCigarType = cigarType;
            currentCigarLength = 1;
        }

        if (base1 != '-')
            ++readBases;
        if (base2 != '-')
            ++refBases;
    }
    alignmentEndPos = alignmentLength;

    m_readEndPos = readBases;
    m_refEndPos = refBases;
    if (currentCigarType == INSERTION && !goToEndSeq1) {
        currentCigarType = CLIP;
        m_readEndPos -= currentCigarLength;
        alignmentEndPos -= currentCigarLength;
    }
    else if (currentCigarType == DELETION && !goToEndSeq2) {
        currentCigarType = NOTHING;
        m_refEndPos -= currentCigarLength;
        alignmentEndPos -= currentCigarLength;
    }

    cigarTypes.push_back(currentCigarType);
    cigarLengths.push_back(currentCigarLength);

    // Build the CIGAR string and tally up the score.
    m_cigar = "";
    int alignmentPos = 0;
    for (size_t i = 0; i < cigarTypes.size(); ++i) {
        CigarType type = cigarTypes[i];
        int length = cigarLengths[i];

        std::string cigarPart = getCigarPart(type, length);
        m_cigar += cigarPart;
        int score = getCigarScore(type, length, scoringScheme, readAlignment, refAlignment, alignmentPos);
        m_rawScore += score;
        alignmentPos += length;
    }
    int alignmentLengthExcludingClips = alignmentEndPos - alignmentStartPos;
    int perfectScore = scoreMatch(scoringScheme) * alignmentLengthExcludingClips;
    int worstScore = scoreMismatch(scoringScheme) * alignmentLengthExcludingClips;
    if (perfectScore > worstScore)
        m_scaledScore = 100.0 * double(m_rawScore - worstScore) / double(perfectScore - worstScore);
    else
        m_scaledScore = 0.0;

    // Add the offset to the reference positions so they reflect the actual alignment location in
    // the full reference, not the trimmed reference that was given to Seqan.
    m_refStartPos += refOffset;
    m_refEndPos += refOffset;

    m_milliseconds = getTime() - startTime;
}


std::string ScoredAlignment::getFullString() {
    std::string revCompStr;
    if (isRevComp())
        revCompStr = '-';
    else
        revCompStr = '+';

    return m_refName + "," + 
           revCompStr + "," + 
           std::to_string(m_readStartPos) + "," + 
           std::to_string(m_readEndPos) + "," + 
           std::to_string(m_refStartPos) + "," + 
           std::to_string(m_refEndPos) + "," + 
           std::to_string(m_rawScore) + "," +
           std::to_string(m_scaledScore) + "," +
           std::to_string(m_milliseconds) + "," +
           m_cigar;
}


std::string ScoredAlignment::getShortDisplayString() {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << m_scaledScore;

    return m_readName + " (" + std::to_string(m_readStartPos) + "-" +  std::to_string(m_readEndPos) + "), " +
           m_refName + " (" + std::to_string(m_refStartPos) + "-" +  std::to_string(m_refEndPos) + "), " +
           "raw score = " + std::to_string(m_rawScore) + ", scaled score = " + ss.str() + ", band size = " + std::to_string(m_bandSize);
}



CigarType ScoredAlignment::getCigarType(char b1, char b2, bool alignmentStarted) {
    if (b1 == '-') {
        if (alignmentStarted)
            return DELETION;
        else
            return NOTHING;
    }
    else if (b2 == '-') {
        if (alignmentStarted)
            return INSERTION;
        else
            return CLIP;
    }
    else
        return MATCH;
}

std::string ScoredAlignment::getCigarPart(CigarType type, int length) {
    std::string cigarPart = std::to_string(length);
    if (type == DELETION)
        cigarPart.append("D");
    else if (type == INSERTION)
        cigarPart.append("I");
    else if (type == CLIP)
        cigarPart.append("S");
    else if (type == MATCH)
        cigarPart.append("M");
    else //type == NOTHING
        return "";
    return cigarPart;
}


int ScoredAlignment::getCigarScore(CigarType type, int length, Score<int, Simple> & scoringScheme,
                                       std::string & readAlignment, std::string & refAlignment,
                                       int alignmentPos) {

    // Scoring indels is easy because we only need to know the length.
    if (type == INSERTION || type == DELETION)
        return scoreGapOpen(scoringScheme) + ((length - 1) * scoreGapExtend(scoringScheme));

    // To score matches we must actually look at the bases.
    else if (type == MATCH) {
        int score = 0;
        for (int i = 0; i < length; ++i) {
            int pos = alignmentPos + i;
            bool match = (readAlignment[pos] == refAlignment[pos]);
            if (match)
                score += scoreMatch(scoringScheme);
            else
                score += scoreMismatch(scoringScheme);
        }
        return score;
    }
    return 0;
}

bool ScoredAlignment::isRevComp() {
    return m_readName.back() == '-';
}


long long getTime() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}