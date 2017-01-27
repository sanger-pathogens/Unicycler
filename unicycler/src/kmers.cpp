// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "kmers.h"


CommonKmer::CommonKmer(int hPosition, int vPosition) :
    m_hPosition(hPosition),
    m_vPosition(vPosition)
{
}


// This is the destructor for KmerPositions. It cleans up all the KmerPosMaps which were allocated
// on the heap.
KmerPositions::~KmerPositions() {
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        delete i->second;
}

// Returns a vector all of k-mer position names (should be exact the same as the the sequence names).
std::vector<std::string> KmerPositions::getAllNames() {
    std::vector<std::string> returnVector;
    m_mutex.lock();
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        returnVector.push_back(i->first);
    m_mutex.unlock();
    return returnVector;
}

// Returns the length of the sequence with the given name.
int KmerPositions::getLength(std::string & name) {
    int returnedLength = 0;
    m_mutex.lock();
    if (m_sequences.find(name) != m_sequences.end())
        returnedLength = m_sequences[name].length();
    m_mutex.unlock();
    return returnedLength;
}

// This function adds a sequence to the KmerPositions object. It creates a new KmerPosMap on the
// heap (will be deleted in destructor), fills it up and adds it to m_kmerPositions.
void KmerPositions::addPositions(std::string & name, std::string & sequence, int kSize) {
    KmerPosMap * posMap = new KmerPosMap();
    int kCount = sequence.size() - kSize + 1;
    for (int i = 0; i < kCount; ++i) {
        std::string kmer = sequence.substr(i, kSize);
        if (posMap->find(kmer) == posMap->end())
            (*posMap)[kmer] = std::vector<int>();
        (*posMap)[kmer].push_back(i);
    }

    m_mutex.lock();
    m_sequences[name] = sequence;
    m_kmerPositions[name] = posMap;
    m_mutex.unlock();
}

// This function retrieves a KmerPosMap from the object using the name as a key. If the name isn't
// in the map, it returns 0.
KmerPosMap * KmerPositions::getKmerPositions(std::string & name) {
    KmerPosMap * returnedMap = 0;
    m_mutex.lock();
    if (m_kmerPositions.find(name) != m_kmerPositions.end())
        returnedMap =  m_kmerPositions[name];
    m_mutex.unlock();
    return returnedMap;
}

std::string * KmerPositions::getSequence(std::string & name) {
    std::string * returnedSequence = 0;
    m_mutex.lock();
    if (m_kmerPositions.find(name) != m_kmerPositions.end())
        returnedSequence = &(m_sequences[name]);
    m_mutex.unlock();
    return returnedSequence;
}

KmerPositions * newKmerPositions() {
    return new KmerPositions();
}

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC, int kSize) {
    std::string name(nameC);
    std::string sequence(sequenceC);
    kmerPositions->addPositions(name, sequence, kSize);
}

void deleteAllKmerPositions(KmerPositions * kmerPositions) {
    delete kmerPositions;
}
