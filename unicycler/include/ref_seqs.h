// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef REF_SEQS_H
#define REF_SEQS_H

#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::string> SeqMap;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    SeqMap * newRefSeqs();
    void addRefSeq(SeqMap * seqMap, char * nameC, char * sequenceC);
    void deleteRefSeqs(SeqMap * seqMap);
}

#endif // REF_SEQS_H
