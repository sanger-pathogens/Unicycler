// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#ifndef SEMI_GLOBAL_ALIGN_H
#define SEMI_GLOBAL_ALIGN_H

#include <seqan/sequence.h>
#include <seqan/seeds.h>
#include <seqan/align.h>
#include <seqan/basic.h>
#include <string>
#include <vector>
#include <unordered_set>
#include "kmers.h"
#include "scoredalignment.h"
#include "random_alignments.h"
#include "string_functions.h"
#include "ref_seqs.h"
#include "nanoflann.hpp"

using namespace seqan;
using namespace nanoflann;

typedef std::pair<int, int> StartEndRange;
typedef std::unordered_map<std::string, std::vector<StartEndRange> > RefRangeMap;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;

struct Point {
    int x, y;
    Point() {x = 0; y = 0;}
    Point(int p_x, int p_y) {x = p_x; y = p_y;}
    bool operator==(const Point &other) const {return x == other.x && y == other.y;}
    bool operator<(const Point &other) const {  // sort using x first, then y
        if (x == other.x)
            return y < other.y;
        return x < other.x;
    }
};

typedef std::unordered_set<Point> PointSet;
typedef std::vector<Point> PointVector;


// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                               char * minimapAlignmentsStr, SeqMap * refSeqs,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore, double lowScoreThreshold, bool returnBad,
                               int sensitivityLevel);
}

std::vector<ScoredAlignment *> alignReadToReferenceRange(SeqMap * refSeqs, std::string refName,
                                                         StartEndRange refRange, int refLen,
                                                         std::string readName, char readStrand,
                                                         KmerPosMap * kmerPositions, int kSize,
                                                         std::string * readSeq, int matchScore,
                                                         int mismatchScore, int gapOpenScore,
                                                         int gapExtensionScore,
                                                         int sensitivityLevel,
                                                         int verbosity, std::string & output);

std::pair<int,int> getRefRange(int refStart, int refEnd, int refLen,
                               int readStart, int readEnd, int readLen, bool posStrand);

std::vector<std::pair<int, int> > simplifyRanges(std::vector<std::pair<int, int> > & ranges);

// http://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key#
namespace std {
  template <>
  struct hash<Point> {
    size_t operator()(const Point& p) const {
      return (std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1)) >> 1;
    }
  };
}

// This stuff is for the nanoflann NN searching.
struct PointCloud
{
    PointVector pts;
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline int kdtree_distance(const int *p1, const size_t idx_p2,size_t /*size*/) const
    {
        const int d0=p1[0]-pts[idx_p2].x;
        const int d1=p1[1]-pts[idx_p2].y;
        return d0+d1;
    }
    inline int kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0)
            return pts[idx].x;
        else
            return pts[idx].y;
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

typedef KDTreeSingleIndexAdaptor<L1_Adaptor<int, PointCloud>, PointCloud, 2> my_kd_tree_t;

PointVector radiusSearchAroundPoint(Point point, int radius, PointCloud & cloud,
                                    my_kd_tree_t & index);

Point getHighestDensityPoint(int densityRadius, PointCloud & cloud, my_kd_tree_t & index,
                             std::string & trimmedRefSeq, std::string * readSeq);

double getPointDensityScore(int densityRadius, Point p, PointCloud & cloud, my_kd_tree_t & index);

void addKmerPointsToNanoflann(PointCloud & cloud, std::vector<CommonKmer> & commonKmers,
                              PointSet & usedPoints);

double getSlope(Point & p1, Point & p2);

PointSet lineTracingWithNanoflann(std::vector<CommonKmer> & commonKmers, PointSet & usedPoints,
                                  PointCloud & cloud, my_kd_tree_t & index, std::string readName,
                                  char readStrand, std::string * readSeq, int readLen,
                                  std::string refName, std::string & trimmedRefSeq, int lineNum,
                                  int verbosity, std::string & output, bool & failedLine,
                                  double & pointSetScore);

void displayRFunctions(std::string & output);

void displayRefRanges(std::string & output, RefRangeMap & simplifiedRefRanges);

void saveCommonKmersToFile(std::string readName, char readStrand, std::string refName,
                           std::vector<CommonKmer> & commonKmers, std::string & output);

void saveChainedSeedsToFile(std::string readName, char readStrand, std::string refName,
                            String<TSeed> & seedChain, std::string & output, int maxLineNum,
                            int bestLineNum);

void saveTraceDotsToFile(std::string readName, char readStrand, std::string refName,
                         PointVector & traceDots, PointSet & pointSet, std::string & output,
                         int lineNum);

double scorePointSet(PointSet & pointSet, PointVector & traceDots, bool & failedLine);

double getWorstSlope(PointVector traceDots);

Point mutateLineToBestFitPoints(Point previousP, Point newP, PointCloud & cloud,
                                my_kd_tree_t & index, PointSet & pointsNearLine,
                                bool leftAlignmentRectangle);

void addPointsNearLine(Point p1, Point p2, PointSet & pointsNearLine, PointSet & pointSet,
                       double radius);

double distanceToLineSegment(Point p, Point l1, Point l2);

double scoreLineSegment(Point p1, Point p2, PointSet & pointsNearLine);

double variance(std::vector<double> & v);

long long getMaxSeedChainGapArea(String<TSeed> & seedChain, int readLen, int trimmedRefLen);

#endif // SEMI_GLOBAL_ALIGN_H
