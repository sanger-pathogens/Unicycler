// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef ALIGNMENT_H
#define ALIGNMENT_H


#include <string>
#include <seqan/basic.h>
#include <seqan/align.h>

using namespace seqan;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};


class ScoredAlignment {
public:
    ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment, 
                    std::string & readName, std::string & refName,
                    int readLength, int refLength,
                    int refOffset, long long startTime, int bandSize,
                    bool startImmediately, bool goToEndSeq1, bool goToEndSeq2,
                    Score<int, Simple> & scoringScheme);
    std::string getFullString();
    std::string getShortDisplayString();
    bool isRevComp();
    int getReadAlignmentLength() {return m_readEndPos - m_readStartPos;}
    int getRefAlignmentLength() {return m_refEndPos - m_refStartPos;}

    std::string m_readName;
    std::string m_refName;
    int m_readLength;
    int m_refLength;
    int m_readStartPos;
    int m_readEndPos;
    int m_refStartPos;
    int m_refEndPos;
    std::string m_cigar;
    int m_rawScore;
    double m_scaledScore;
    int m_milliseconds;
    int m_bandSize;

private:
    CigarType getCigarType(char b1, char b2, bool alignmentStarted);
    std::string getCigarPart(CigarType type, int length);
    int getCigarScore(CigarType type, int length, Score<int, Simple> & scoringScheme,
                      std::string & readAlignment, std::string & refAlignment,
                      int alignmentPos);
};

long long getTime();

#endif // ALIGNMENT_H
