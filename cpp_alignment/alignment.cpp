

Alignment::Alignment(Align<Dna5String, ArrayGaps> & alignment, int refOffset, long long startTime,
                     bool startImmediately, bool goToEnd, ScoringScheme & scoringScheme):
    m_readStartPos(-1), m_refStartPos(-1), m_rawScore(0) {
    // Extract the alignment sequences into C++ strings, as the TRow type doesn't seem to have
    // constant time random access.
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

    if (startImmediately) {
        alignmentStarted = true;
        readStarted = true;
        refStarted = true;
        m_readStartPos = 0;
        m_refStartPos = 0;
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
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (cigarType == MATCH) {
            if (base1 == base2)
                m_rawScore += scoreMatch(scoringScheme);
            else
                m_rawScore += scoreMismatch(scoringScheme);
        }
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

    m_readEndPos = readBases;
    m_refEndPos = refBases;
    if (currentCigarType == INSERTION && !goToEnd) {
        currentCigarType = CLIP;
        m_readEndPos -= currentCigarLength;
    }
    else if (currentCigarType == DELETION && !goToEnd) {
        currentCigarType = NOTHING;
        m_refEndPos -= currentCigarLength;
    }

    cigarTypes.push_back(currentCigarType);
    cigarLengths.push_back(currentCigarLength);

    // Build the CIGAR string and tally up indel scores.
    m_cigar = "";
    for (int i = 0; i < cigarTypes.size(); ++i) {
        m_cigar += getCigarPart(cigarTypes[i], cigarLengths[i]);
        m_rawScore += getCigarScore(cigarTypes[i], cigarLengths[i]);
    }
    int perfectScore = scoreMatch(scoringScheme) * alignmentLength;
    m_scaledScore = 100.0 * double(m_rawScore) / perfectScore;

    m_milliseconds = getTime() - startTime;
}


std::string Alignment::getFullString() {
    return m_cigar + ";" +
           std::to_string(m_readStartPos) + "," + 
           std::to_string(m_readEndPos) + "," + 
           std::to_string(m_refStartPos) + "," + 
           std::to_string(m_refEndPos) + "," + 
           std::to_string(m_milliseconds);
}


std::string Alignment::getShortDisplayString() {
    return "(" + std::to_string(m_readStartPos) + "-" +  std::to_string(m_readEndPos) + "), " + 
           "(" + std::to_string(m_refStartPos) + "-" +  std::to_string(m_refEndPos) + "), " + 
           "raw score = " + std::to_string(m_rawScore) + ", scaled score = " + std::to_string(m_scaledScore)
}



CigarType Alignment::getCigarType(char b1, char b2, bool alignmentStarted) {
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

std::string Alignment::getCigarPart(CigarType type, int length) {
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

// Only returns scores for indel cigar parts, as match/mismatch scores are done on a base-by-base
// basis.
int getCigarScore(CigarType type, int length, ScoringScheme & scoringScheme) {
    if (type == INSERTION || type == DELETION)
        return scoreGapOpen(scoringScheme) + ((length - 1) * scoreGapExtend(scoringScheme);
    else
        return 0;
}



int Alignment::getIndelScore(int length, ScoringScheme & scoringScheme) {
    return scoreGapOpen(scoringScheme) + ((currentCigarLength - 1) * scoreGapExtend(scoringScheme));
}
