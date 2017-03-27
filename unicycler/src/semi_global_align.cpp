// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#include "semi_global_align.h"

#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <utility>
#include <math.h>

#include "settings.h"


char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                           char * minimapAlignmentsStr, SeqMap * refSeqs,
                           int matchScore, int mismatchScore, int gapOpenScore,
                           int gapExtensionScore, double /*lowScoreThreshold*/, bool /*returnBad*/,
                           int sensitivityLevel) {
    int kSize = LEVEL_0_KMER_SIZE;
    if (sensitivityLevel == 1)
        kSize = LEVEL_1_KMER_SIZE;
    else if (sensitivityLevel == 2)
        kSize = LEVEL_2_KMER_SIZE;
    else if (sensitivityLevel == 3)
        kSize = LEVEL_3_KMER_SIZE;

    std::string output;
    std::string returnString;
    std::vector<ScoredAlignment *> returnedAlignments;

    // Change the read name and sequence to C++ strings.
    std::string readName(readNameC);
    std::string posReadName = readName + "+";
    std::string negReadName = readName + "-";
    std::string posReadSeq(readSeqC);
    std::string negReadSeq;  // Will make later, if necessary.
    int readLength = int(posReadSeq.length());

    std::vector<std::string> minimapAlignments = splitString(minimapAlignmentsStr, ';');
    if (verbosity > 2) {
        output += "minimap alignments:\n";
        for (auto const & minimapAlignment : minimapAlignments)
            output += "    " + minimapAlignment + "\n";
    }
    if (verbosity > 3)
        displayRFunctions(output);

    // For each minimap alignment we find the appropriate part of the reference sequence.
    RefRangeMap refRanges;
    for (size_t i = 0; i < minimapAlignments.size(); ++i) {
        std::string minimapStr = minimapAlignments[i];
        std::vector<std::string> minimapStrParts = splitString(minimapStr, ',');

        int readStart = std::stoi(minimapStrParts[0]);
        int readEnd = std::stoi(minimapStrParts[1]);
        char readStrand = minimapStrParts[2][0];
        bool posStrand = readStrand == '+';

        std::string refName = minimapStrParts[3];
        int refStart = std::stoi(minimapStrParts[4]);
        int refEnd = std::stoi(minimapStrParts[5]);
        std::string & refSeq = refSeqs->at(refName);
        int refLength = int(refSeq.length());

        StartEndRange refRange = getRefRange(refStart, refEnd, refLength, readStart, readEnd,
                                             readLength, posStrand);

        // The first time we see the reference/strand, initialise it with an empty vector.
        std::string refNameAndStrand = refName + readStrand;
        if (refRanges.find(refNameAndStrand) == refRanges.end())
            refRanges[refNameAndStrand] = std::vector<StartEndRange>();

        refRanges[refNameAndStrand].push_back(refRange);
    }

    // Simplify the reference ranges by combining overlapping ranges.
    RefRangeMap simplifiedRefRanges;
    for(auto const & r : refRanges) {
        std::string refNameAndStrand = r.first;
        std::vector<StartEndRange> ranges = r.second;
        std::vector<StartEndRange> simplifiedRanges = simplifyRanges(ranges);
        simplifiedRefRanges[refNameAndStrand] = simplifiedRanges;
    }
    if (verbosity > 2)
        displayRefRanges(output, simplifiedRefRanges);

    // Make a new KmerPositions object for the read. We'll actually add positions later as
    // necessary (because we may not need both the positive strand or the negative strand).
    KmerPositions readKmerPositions;
    bool posPositions = false, negPositions = false;

    // Align to each reference range.
    for(auto const & r : simplifiedRefRanges) {
        std::string refName = r.first;
        char readStrand = refName.back();
        bool posStrand = readStrand == '+';
        refName.pop_back();
        std::string & refSeq = refSeqs->at(refName);
        int refLength = int(refSeq.length());
        std::vector<StartEndRange> ranges = r.second;

        // Prepare some stuff for the read.
        std::string * readSeq;
        KmerPosMap * kmerPositions;
        if (posStrand) {
            if (!posPositions) {
                readKmerPositions.addPositions(posReadName, posReadSeq, kSize);
                posPositions = true;
            }
            readSeq = &posReadSeq;
            kmerPositions = readKmerPositions.getKmerPositions(posReadName);
        }
        else {  // negative strand
            if (!negPositions) {
                negReadSeq = getReverseComplement(posReadSeq);
                readKmerPositions.addPositions(negReadName, negReadSeq, kSize);
                negPositions = true;
            }
            readSeq = &negReadSeq;
            kmerPositions = readKmerPositions.getKmerPositions(negReadName);
        }

        // Work on each range (there's probably just one, but there could be more).
        for (auto const & range : ranges) {
            std::vector<ScoredAlignment *> a =
                alignReadToReferenceRange(refSeqs, refName, range, refLength, readName, readStrand,
                                          kmerPositions, kSize, readSeq, matchScore, mismatchScore,
                                          gapOpenScore, gapExtensionScore, sensitivityLevel,
                                          verbosity, output);
            returnedAlignments.insert(returnedAlignments.end(), a.begin(), a.end());
        }
    }

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    for (auto const & alignment : returnedAlignments) {
        if (alignment != 0)
            returnString += alignment->getFullString() + ";";
    }
    returnString += output;

    return cppStringToCString(returnString);
}


std::vector<ScoredAlignment *> alignReadToReferenceRange(SeqMap * refSeqs, std::string refName,
                                                         StartEndRange refRange, int refLen,
                                                         std::string readName, char readStrand,
                                                         KmerPosMap * kmerPositions, int kSize,
                                                         std::string * readSeq, int matchScore,
                                                         int mismatchScore, int gapOpenScore,
                                                         int gapExtensionScore,
                                                         int sensitivityLevel,
                                                         int verbosity, std::string & output) {
    long long startTime = getTime();

    // Set parameters based on the sensitivity level.
    int bandSize = LEVEL_0_BAND_SIZE;
    int minLineTraceCount = LEVEL_0_MIN_LINE_TRACE_COUNT;
    int maxLineTraceCount = LEVEL_0_MAX_LINE_TRACE_COUNT;
    if (sensitivityLevel == 1) {
        bandSize = LEVEL_1_BAND_SIZE;
        minLineTraceCount = LEVEL_1_MIN_LINE_TRACE_COUNT;
        maxLineTraceCount = LEVEL_1_MAX_LINE_TRACE_COUNT;
    }
    else if (sensitivityLevel == 2) {
        bandSize = LEVEL_2_BAND_SIZE;
        minLineTraceCount = LEVEL_2_MIN_LINE_TRACE_COUNT;
        maxLineTraceCount = LEVEL_2_MAX_LINE_TRACE_COUNT;
    }
    else if (sensitivityLevel == 3) {
        bandSize = LEVEL_3_BAND_SIZE;
        minLineTraceCount = LEVEL_3_MIN_LINE_TRACE_COUNT;
        maxLineTraceCount = LEVEL_3_MAX_LINE_TRACE_COUNT;
    }

    // Extract the part of the reference to which we're aligning the read.
    int refStart = refRange.first;
    int refEnd = refRange.second;
    int readLen = int(readSeq->length());
    std::string trimmedRefSeq = refSeqs->at(refName).substr(size_t(refStart),
                                                            size_t(refEnd - refStart));
    int trimmedRefLen = int(trimmedRefSeq.length());
    if (verbosity > 2)
        output += "Range: " + refName + ": " + std::to_string(refStart) + " - " + std::to_string(refEnd) + "\n";

    // Find all common k-mer positions.
    std::vector<CommonKmer> commonKmers;
    int maxI = trimmedRefLen - kSize + 1;
    for (int i = 0; i < maxI; ++i) {
        std::string refKmer = trimmedRefSeq.substr(size_t(i), size_t(kSize));
        if (kmerPositions->find(refKmer) != kmerPositions->end() ) {  // if k-mer is in the read
            std::vector<int> & readPositions = kmerPositions->at(refKmer);
            for (size_t j = 0; j < readPositions.size(); ++j)
                commonKmers.emplace_back(readPositions[j], i);
        }
    }
    if (verbosity > 2)
        output += "    common " + std::to_string(kSize) + "-mers: " + std::to_string(commonKmers.size()) + "\n";
    if (verbosity > 3)
        saveCommonKmersToFile(readName, readStrand, refName, commonKmers, output);

    // Build a nanoflann point cloud with all of the common k-mer points.
    PointSet usedPoints;
    PointCloud cloud;
    addKmerPointsToNanoflann(cloud, commonKmers, usedPoints);
    my_kd_tree_t index(2, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    // Use nanoflann and line tracing to get a set of common k-mer positions around a line.
    PointSet bestPointSet;
    double bestPointScore = 0.0;
    int bestLineNum = 0;
    int maxLineNum = 0;
    for (int lineNum = 0; lineNum < maxLineTraceCount; ++lineNum) {
        maxLineNum = lineNum;
        bool failedLine = false;
        double pointSetScore = 0.0;
        PointSet pointSet = lineTracingWithNanoflann(commonKmers, usedPoints, cloud, index,
                                                     readName, readStrand, readSeq, readLen,
                                                     refName, trimmedRefSeq, lineNum, verbosity,
                                                     output, failedLine, pointSetScore);
        if (pointSetScore > bestPointScore) {
            bestPointSet = pointSet;
            bestPointScore = pointSetScore;
            bestLineNum = lineNum;
        }

        // If this line looked good (i.e. we didn't get 'lost' in the line tracing) and we've tried
        // enough lines, then we're done!
        if (!failedLine && lineNum >= minLineTraceCount - 1)
            break;

        // If we've used all the points, we can't do another line!
        if (usedPoints.size() >= commonKmers.size())
            break;
    }

    // Now add the points to Seqan and get a global chain so we can do a banded alignment. We sort
    // them first so they're added in a consistent order.
    String<TSeed> seeds;
    PointVector bestPointSetVector;
    bestPointSetVector.reserve(bestPointSet.size());
    for (auto const & p : bestPointSet)
        bestPointSetVector.push_back(p);
    std::sort(bestPointSetVector.begin(), bestPointSetVector.end());
    for (auto const & p : bestPointSetVector)
        appendValue(seeds, TSeed(size_t(p.x), size_t(p.y), size_t(kSize)));
    TSeedSet seedSet;
    for (unsigned i = 0; i < length(seeds); ++i) {
        if (!addSeed(seedSet, seeds[i], 2, Merge()))
            addSeed(seedSet, seeds[i], Single());
    }
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    if (verbosity > 3)
        saveChainedSeedsToFile(readName, readStrand, refName, seedChain, output, maxLineNum,
                               bestLineNum);

    // If the seed chain contains too much gap area, then we don't proceed - it would take too long
    // to align and is probably not a good alignment anyway.
    std::vector<ScoredAlignment *> alignments;
    int seedChainLength = length(seedChain);
    if (seedChainLength == 0)
        return alignments;
    long long gapArea = getMaxSeedChainGapArea(seedChain, readLen, trimmedRefLen);
    if (gapArea > MAX_BANDED_ALIGNMENT_GAP_AREA)
        return alignments;

    // Finally we can actually do the Seqan alignment!
    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), *readSeq);
    assignSource(row(alignment, 1), trimmedRefSeq);
    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
    ScoredAlignment * sgAlignment;
    try {
        bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig,
                             (unsigned int)bandSize);
        std::string signedReadName = readName + readStrand;
        sgAlignment = new ScoredAlignment(alignment, signedReadName, refName, readLen, refLen,
                                          refStart, startTime, bandSize, false, false, false,
                                          scoringScheme);
        alignments.push_back(sgAlignment);
    }
    catch (...) {}

    return alignments;
}


// This function takes a Seqan seed chain and returns the area of the largest gap. This is because
// Seqan will do full alignments in those gaps, and so big gaps will be slow and take lots of
// memory.
long long getMaxSeedChainGapArea(String<TSeed> & seedChain, int readLen, int trimmedRefLen) {
    int seedChainLength = length(seedChain);
    int previousH = 0;
    int previousV = 0;
    long long maxGapArea = 0;
    for (int i = 0; i <= seedChainLength; ++i) {
        int hPos, vPos;
        if (i == seedChainLength) {
            hPos = readLen;
            vPos = trimmedRefLen;
        }
        else {
            hPos = beginPositionH(seedChain[i]);
            vPos = beginPositionV(seedChain[i]);
        }
        long long hGap = hPos - previousH;
        long long vGap = vPos - previousV;
        long long gapArea = hGap * vGap;
        if (gapArea > maxGapArea)
            maxGapArea = gapArea;
        previousH = hPos;
        previousV = vPos;
    }
    return maxGapArea;
}


PointSet lineTracingWithNanoflann(std::vector<CommonKmer> & commonKmers, PointSet & usedPoints,
                                  PointCloud & cloud, my_kd_tree_t & index, std::string readName,
                                  char readStrand, std::string * readSeq, int readLen,
                                  std::string refName, std::string & trimmedRefSeq, int lineNum,
                                  int verbosity, std::string & output, bool & failedLine,
                                  double & pointSetScore) {
    // First find the highest density point in the region, which we will use to start the trace.
    PointCloud startingPointCloud;
    addKmerPointsToNanoflann(startingPointCloud, commonKmers, usedPoints);
    my_kd_tree_t startingPointIndex(2, startingPointCloud, KDTreeSingleIndexAdaptorParams(10));
    startingPointIndex.buildIndex();
    Point startPoint = getHighestDensityPoint(LINE_TRACING_START_POINT_SEARCH_RADIUS,
                                              startingPointCloud, startingPointIndex,
                                              trimmedRefSeq, readSeq);
    Point p = startPoint;
    PointVector traceDots;
    traceDots.push_back(p);

    // Set the points around the starting point as 'used' so they are excluded from future starting
    // points (to ensure that subsequent line traces begin from a sufficiently different location).
    PointVector nearStart = radiusSearchAroundPoint(p, 2 * LINE_TRACING_START_POINT_SEARCH_RADIUS,
                                                        cloud, index);
    usedPoints.insert(nearStart.begin(), nearStart.end());

    // Start the point collection using points around the starting point.
    PointVector nearbyPoints = radiusSearchAroundPoint(p, TRACE_LINE_COLLECTION_DISTANCE, cloud,
                                                       index);
    PointSet pointSet(nearbyPoints.begin(), nearbyPoints.end());

    // Trace the line forward then backward.
    int directions[2] = {1, -1};
    for (auto const & direction : directions) {
        p = startPoint;
        int maxX = readLen;
        int maxY = int(trimmedRefSeq.length());
        while (true) {
            int step = direction * TRACE_LINE_STEP_DISTANCE;
            Point previousP = p;
            Point newP(p.x + step, p.y + step);

            PointSet pointsNearLine;
            p = mutateLineToBestFitPoints(previousP, newP, cloud, index, pointsNearLine);
            traceDots.push_back(p);
            addPointsNearLine(previousP, p, pointsNearLine, pointSet,
                              TRACE_LINE_COLLECTION_DISTANCE);

            // Quit advancing when we've left the alignment rectangle.
            if (direction == 1 && (p.x > maxX || p.y > maxY))
                break;
            if (direction == -1 && (p.x < 0 || p.y < 0))
                break;
        }
    }
    pointSetScore = scorePointSet(pointSet, traceDots, failedLine);

    if (verbosity > 2) {
        output += "    line " + std::to_string(lineNum + 1) + ": ";
        output += std::to_string(pointSet.size()) + " points, ";
        output += "score=" + std::to_string(pointSetScore) + " (";
        if (failedLine)
            output += "bad";
        else
            output += "good";
        output += ")\n";
    }
    if (verbosity > 3)
        saveTraceDotsToFile(readName, readStrand, refName, traceDots, pointSet, output, lineNum);
    return pointSet;
}


Point shiftPointUp(Point p, int steps) {
    p.x -= steps;
    p.y += steps;
    return p;
}


Point shiftPointDown(Point p, int steps) {
    p.x += steps;
    p.y -= steps;
    return p;
}


Point mutateLineToBestFitPoints(Point p1, Point p2, PointCloud & cloud, my_kd_tree_t & index,
                                PointSet & pointsNearLine) {

    int radius = int(TRACE_LINE_STEP_DISTANCE * 1.1);
    PointVector pointsNearP1 = radiusSearchAroundPoint(p1, radius, cloud, index);
    PointVector pointsNearP2 = radiusSearchAroundPoint(p2, radius, cloud, index);
    pointsNearLine.insert(pointsNearP1.begin(), pointsNearP1.end());
    pointsNearLine.insert(pointsNearP2.begin(), pointsNearP2.end());

    Point p2Up = shiftPointUp(p2, TRACE_LINE_MUTATION_SIZE);
    Point p2Down = shiftPointDown(p2, TRACE_LINE_MUTATION_SIZE);
    double unmutatedScore = scoreLineSegment(p1, p2, pointsNearLine);
    double mutatedUpScore = scoreLineSegment(p1, p2Up, pointsNearLine);
    double mutatedDownScore = scoreLineSegment(p1, p2Down, pointsNearLine);

    while (true) {
        // If neither mutation helps, then we're done!
        if (unmutatedScore >= mutatedUpScore && unmutatedScore >= mutatedDownScore)
            break;

        else if (mutatedUpScore > unmutatedScore) {
            p2Down = p2;
            mutatedDownScore = unmutatedScore;
            p2 = p2Up;
            unmutatedScore = mutatedUpScore;
            p2Up = shiftPointUp(p2, TRACE_LINE_MUTATION_SIZE);
            mutatedUpScore = scoreLineSegment(p1, p2Up, pointsNearLine);
        }

        else if (mutatedDownScore > unmutatedScore) {
            p2Up = p2;
            mutatedUpScore = unmutatedScore;
            p2 = p2Down;
            unmutatedScore = mutatedDownScore;
            p2Down = shiftPointDown(p2, TRACE_LINE_MUTATION_SIZE);
            mutatedDownScore = scoreLineSegment(p1, p2Down, pointsNearLine);
        }
    }
    return p2;
}


// Line segments are scored on two fronts: their slope (closer to 1 is better) and the closeness
// of points to the line segment.
double scoreLineSegment(Point p1, Point p2, PointSet & pointsNearLine) {
    double slope = getSlope(p1, p2);
    if (slope > 1.0)
        slope = 1.0 / slope;
    double slopeScore = (MAX_SLOPE_SCORE / (1.0 - MIN_ACCEPTABLE_LINE_SEGMENT_SLOPE)) *
                        (slope - MIN_ACCEPTABLE_LINE_SEGMENT_SLOPE);
    double maxScorePerPoint = MAX_POINTS_SCORE / TRACE_LINE_STEP_DISTANCE;
    double pointDistanceScore = 0.0;
    for (auto const & p : pointsNearLine) {
        double dist = distanceToLineSegment(p, p1, p2);
        pointDistanceScore += maxScorePerPoint / (dist + 1.0);
    }
    double finalScore = slopeScore + pointDistanceScore;
    return finalScore;
}

void addPointsNearLine(Point p1, Point p2, PointSet & pointsNearLine, PointSet & pointSet,
                       double radius) {
    for (auto const & p : pointsNearLine) {
        if (distanceToLineSegment(p, p1, p2) <= radius)
            pointSet.insert(p);
    }
}


// http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
double distanceToLineSegment(Point p, Point l1, Point l2) {
    double A = p.x - l1.x;
    double B = p.y - l1.y;
    double C = l2.x - l1.x;
    double D = l2.y - l1.y;
    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1;
    if (len_sq != 0) //in case of 0 length line
        param = dot / len_sq;
    double xx, yy;
    if (param < 0) {
        xx = l1.x;
        yy = l1.y;
    }
    else if (param > 1) {
        xx = l2.x;
        yy = l2.y;
    }
    else {
        xx = l1.x + param * C;
        yy = l1.y + param * D;
    }
    double dx = p.x - xx;
    double dy = p.y - yy;
    return sqrt(dx * dx + dy * dy);
}


void addKmerPointsToNanoflann(PointCloud & cloud, std::vector<CommonKmer> & commonKmers,
                              PointSet & usedPoints) {
    PointVector points;
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        Point p(commonKmers[i].m_hPosition, commonKmers[i].m_vPosition);
        bool alreadyUsed = usedPoints.find(p) != usedPoints.end();
        if (!alreadyUsed)
            points.push_back(p);
    }
    size_t pCount = points.size();
    cloud.pts.resize(pCount);
    for (size_t i = 0; i < pCount; ++i) {
        cloud.pts[i].x = points[i].x;
        cloud.pts[i].y = points[i].y;
    }
}


PointVector radiusSearchAroundPoint(Point point, int radius, PointCloud & cloud,
                                    my_kd_tree_t & index) {
    PointVector points;
    nanoflann::SearchParams params;
    std::vector<std::pair<size_t,int> > ret_matches;
    const int query_pt[2] = {point.x, point.y};
    index.radiusSearch(query_pt, radius, ret_matches, params);
    for (auto const & i : ret_matches)
        points.push_back(cloud.pts[i.first]);
    return points;
}



Point getHighestDensityPoint(int densityRadius, PointCloud & cloud, my_kd_tree_t & index,
                             std::string & trimmedRefSeq, std::string * readSeq) {
    PointVector points = getPointsInHighestDensityRegion(densityRadius * 2, trimmedRefSeq, readSeq,
                                                         cloud, index);
    Point highestDensityPoint = points[0];
    double highestDensityScore = 0.0;
    for (auto const & point : points) {
        double densityScore = getPointDensityScore(densityRadius, point, cloud, index);
        if (densityScore > highestDensityScore) {
            highestDensityScore = densityScore;
            highestDensityPoint = point;
        }
    }
    return highestDensityPoint;
}


PointVector getPointsInHighestDensityRegion(int searchRadius, std::string & trimmedRefSeq,
                                            std::string * readSeq, PointCloud & cloud,
                                            my_kd_tree_t & index) {
    int xStepCount = int(ceil(readSeq->length() / double(searchRadius)));
    int yStepCount = int(ceil(trimmedRefSeq.length() / double(searchRadius)));
    double xStepSize = double(readSeq->length()) / xStepCount;
    double yStepSize = double(trimmedRefSeq.length()) / yStepCount;

    nanoflann::SearchParams params;
    double highestDensity = 0.0;
    PointVector pointsInHighestDensity;

    for (int i = 0; i <= xStepCount; ++i) {
        int xCentre = int(0.5 + i * xStepSize);

        for (int j = 0; j <= yStepCount; ++j) {
            int yCentre = int(0.5 + j * yStepSize);

            const int query_pt[2] = {xCentre, yCentre};

            std::vector<std::pair<size_t,int> > ret_matches;
            const size_t nMatches = index.radiusSearch(query_pt, searchRadius, ret_matches, params);
            double density = double(nMatches);

            // Search regions on the edge will have less density than they should (because the
            // region has less area). This biases the search away from the edges, but that's okay
            // because alignments often get tricky and repetitive near the edges (i.e. the ends of
            // contigs) and so we probably don't want to start our line tracing there.

            if (density > highestDensity) {
                highestDensity = density;
                pointsInHighestDensity.clear();
                for (auto const & k : ret_matches)
                    pointsInHighestDensity.push_back(cloud.pts[k.first]);
            }
        }
    }
    return pointsInHighestDensity;
}


// For a given point, the function scores it based on the density of nearby points. Specifically,
// it rewards points that have lots of neighbours close to the diagonal, but it punishes points
// with too many neighbours away from the diagonal.
double getPointDensityScore(int densityRadius, Point p, PointCloud & cloud, my_kd_tree_t & index) {
    PointVector neighbourPoints = radiusSearchAroundPoint(p, densityRadius, cloud, index);
    double a = 1.0 / SCORE_DISTANCE_FROM_DIAGONAL;
    double densityScore = 0.0;
    for (auto const & neighbourPoint : neighbourPoints) {
        int xDiff = neighbourPoint.x - p.x;
        int yDiff = neighbourPoint.y - p.y;
        densityScore += ((1.0 + a) / (abs(xDiff-yDiff) + 1.0)) - a;
    }
    return densityScore;
}


std::pair<int, int> getRefRange(int refStart, int refEnd, int refLen,
                                int readStart, int readEnd, int readLen, bool posStrand) {
    int halfReadLen = 1 + readLen / 2;
    int readBasesBeforeStart = readStart;
    int readBasesAfterEnd = readLen - readEnd;
    if (!posStrand)
        std::swap(readBasesBeforeStart, readBasesAfterEnd);
    int newRefStart = refStart - readBasesBeforeStart - halfReadLen;
    int newRefEnd = refEnd + readBasesAfterEnd + halfReadLen;
    newRefStart = std::max(0, newRefStart);
    newRefEnd = std::min(refLen, newRefEnd);
    return std::pair<int, int>(newRefStart, newRefEnd);
}


std::vector<StartEndRange> simplifyRanges(std::vector<StartEndRange> & ranges) {
    std::sort(ranges.begin(), ranges.end());
    std::vector<std::pair<int, int> > simplifiedRanges;
    std::vector<std::pair<int, int> >::iterator it = ranges.begin();
    std::pair<int,int> current = *(it)++;
    while (it != ranges.end()){
       if (current.second >= it->first){
           current.second = std::max(current.second, it->second); 
       } else {
           simplifiedRanges.push_back(current);
           current = *(it);
       }
       it++;
    }
    simplifiedRanges.push_back(current);
    return simplifiedRanges;
}


double getSlope(Point & p1, Point & p2) {
    int xDiff = p1.x - p2.x;
    int yDiff = p1.y - p2.y;
    double slope = 0.0;
    if (xDiff != 0)
        slope = double(yDiff) / double(xDiff);
    if (xDiff == 0 && yDiff == 0)
        slope = 1.0;
    return slope;
}


void displayRFunctions(std::string & output) {
    output += "R_code:library(ggplot2)\n";
    output += "R_code:library(readr)\n";
    output += "R_code:dot.plot.1 <- function(all_points) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size_1, alpha=p_alpha_1, shape=19) + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
    output += "R_code:dot.plot.2 <- function(all_points, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size_1, alpha=p_alpha_2, shape=19) + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=p_size_2, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
    output += "R_code:dot.plot.3 <- function(all_points, filtered_data, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size_1, alpha=p_alpha_2, shape=19) + geom_point(data=filtered_data,  aes(x=X1, y=X2), size=p_size_1, alpha=1, shape=19, colour=\"green\") + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=p_size_2, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
    output += "R_code:p_size_1 <- 0.1\n";
    output += "R_code:p_size_2 <- 1.0\n";
    output += "R_code:p_alpha_1 <- 0.1\n";
    output += "R_code:p_alpha_2 <- 0.02\n";
}


void displayRefRanges(std::string & output, RefRangeMap & simplifiedRefRanges) {
    output += "Reference ranges:\n";
    for(auto const & r : simplifiedRefRanges) {
        std::string refName = r.first;
        std::vector<StartEndRange> ranges = r.second;
        for (auto const & refRange : ranges) {
            int refStart = refRange.first;
            int refEnd = refRange.second;
            output += "    " + refName + ": " + std::to_string(refStart) + " - " + std::to_string(refEnd) + "\n";
        }
    }
}


void saveCommonKmersToFile(std::string readName, char readStrand, std::string refName,
                           std::vector<CommonKmer> & commonKmers, std::string & output) {
    std::ofstream allPointsFile;
    std::string filename = readName + readStrand + "_" + refName + "_all_points.tsv";
    allPointsFile.open(filename);
    for (auto const & k : commonKmers)
        allPointsFile << k.m_hPosition << "\t" << k.m_vPosition << "\n";
    allPointsFile.close();
    output += "R_code:    all.points <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
}


void saveChainedSeedsToFile(std::string readName, char readStrand, std::string refName,
                            String<TSeed> & seedChain, std::string & output, int maxLineNum,
                            int bestLineNum) {
    std::ofstream chainedSeedsFile;
    std::string filename = readName + readStrand + "_" + refName + "_chained_seeds.tsv";
    chainedSeedsFile.open(filename);
    for (unsigned i = 0; i < length(seedChain); ++i)
        chainedSeedsFile << beginPositionH(seedChain[i]) << "\t" << beginPositionV(seedChain[i]) << "\n";
    chainedSeedsFile.close();
    output += "R_code:    chained.seeds <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
    output += "R_code:    x_limits <- c(0, max(all.points[,1]))\n";
    output += "R_code:    y_limits <- c(0, max(all.points[,2]))\n";
    output += "R_code:    dot.plot.1(all.points)\n";
    for (int lineNum = 0; lineNum <= maxLineNum; ++lineNum) {
        std::string lineNumStr = std::to_string(lineNum+1);
        output += "R_code:    dot.plot.2(all.points, trace.dots." + lineNumStr + ")\n";
        output += "R_code:    dot.plot.3(all.points, filtered.data." + lineNumStr + ", trace.dots." + lineNumStr + ")\n";
    }
    output += "R_code:    dot.plot.3(all.points, chained.seeds, trace.dots." + std::to_string(bestLineNum+1) + ")\n";
}


void saveTraceDotsToFile(std::string readName, char readStrand, std::string refName,
                         PointVector & traceDots, PointSet & pointSet, std::string & output,
                         int lineNum) {
    std::ofstream traceDotsFile;
    std::string lineNumStr = std::to_string(lineNum+1);
    std::string filename = readName + readStrand + "_" + refName + "_line_" + lineNumStr + "_trace_dots.tsv";
    traceDotsFile.open(filename);
    for (auto const & d : traceDots)
        traceDotsFile << d.x << "\t" << d.y << "\n";
    traceDotsFile.close();
    output += "R_code:        trace.dots." + lineNumStr + " <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";

    std::ofstream filteredDataFile;
    filename = readName + readStrand + "_" + refName + "_line_" + lineNumStr + "_filtered_data.tsv";
    filteredDataFile.open(filename);
    for (auto const & d : pointSet)
        filteredDataFile << d.x << "\t" << d.y << "\n";
    filteredDataFile.close();
    output += "R_code:        filtered.data." + lineNumStr + " <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
}



// Standard deviation (http://stackoverflow.com/questions/7616511/)
double variance(std::vector<double> & v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return sq_sum / v.size();
}


// This function gives a point set a quality score so we can choose between alternative point sets
// for an alignment. It also labels the point set as failed or not.
double scorePointSet(PointSet & pointSet, PointVector & traceDots, bool & failedLine) {

    // More points is better.
    double pointCount = double(pointSet.size());

    // Slopes near 1 are better.
    double worstSlope = getWorstSlope(traceDots);
    double worstSlopeScore = worstSlope * 0.9 + 0.1;

    // A good line should have points evenly distributed over its length.
    std::vector<double> xPlusY;
    double min = std::numeric_limits<int>::max();
    double max = std::numeric_limits<int>::min();
    xPlusY.reserve(pointSet.size());
    for (auto const & p : pointSet) {
        double sum = p.x + p.y;
        xPlusY.push_back(sum);
        if (sum > max)
            max = sum;
        if (sum < min)
            min = sum;
    }
    double uniformDistributionVariance = (max - min) * (max - min) / 12.0;
    double thisVariance = variance(xPlusY);
    double varianceScore = thisVariance / uniformDistributionVariance;
    if (varianceScore > 1.0)
        varianceScore = 1.0 / varianceScore;

    // If the slope and variance look week, we say the line has failed.
    failedLine = (worstSlopeScore * varianceScore) < 0.8;

    return pointCount * worstSlopeScore * varianceScore;
}


double getWorstSlope(PointVector traceDots) {
    double worstSlope = 1.0;
    std::sort(traceDots.begin(), traceDots.end());
    for (size_t i = 0; i < traceDots.size() - 1; ++i) {
        Point d1 = traceDots[i];
        Point d2 = traceDots[i + 1];
        double slope = getSlope(d1, d2);
        if (slope > 1)
            slope = 1.0 / slope;
        if (slope < worstSlope)
            worstSlope = slope;
    }
    return worstSlope;
}
