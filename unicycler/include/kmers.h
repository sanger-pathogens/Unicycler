// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef KMERS_H
#define KMERS_H

#include <string>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "settings.h"
#include <mutex>

typedef std::unordered_map<std::string, std::vector<int> > KmerPosMap;


class CommonKmer {
public:
    CommonKmer(int hPosition, int vPosition);
    int m_hPosition;
    int m_vPosition;
};


// KmerPositions is a class that holds maps of k-mer positions for named sequences. It exists so we
// don't have to repeatedly find the same k-mer sets over and over.
class KmerPositions {
public:
    KmerPositions() {}
    ~KmerPositions();
    void addPositions(std::string & name, std::string & sequence, int kSize);
    KmerPosMap * getKmerPositions(std::string & name);
    std::string * getSequence(std::string & name);
    std::vector<std::string> getAllNames();
    int getLength(std::string & name);

private:
    std::unordered_map<std::string, KmerPosMap *> m_kmerPositions;
    std::unordered_map<std::string, std::string> m_sequences;
    std::mutex m_mutex;
};

KmerPositions * newKmerPositions();

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC, int kSize);

void deleteAllKmerPositions(KmerPositions * kmerPositions);


#endif // KMERS_H

