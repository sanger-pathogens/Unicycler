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

#include <seqan/align.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <math.h>

#include "settings.h"
#include <seqan/basic.h>
#include <seqan/seeds.h>


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

//    std::cout << std::endl;  // TEMP

    // Change the read name and sequence to C++ strings.
    std::string readName(readNameC);
    std::string posReadName = readName + "+";
    std::string negReadName = readName + "-";
    std::string posReadSeq(readSeqC);
    std::string negReadSeq;  // Will make later, if necessary.
    int readLength = posReadSeq.length();


    std::vector<std::string> minimapAlignments = splitString(minimapAlignmentsStr, ';');
    if (verbosity > 2) {
        output += "minimap alignments:\n";
        for (auto minimapAlignment : minimapAlignments)
            output += "    " + minimapAlignment + "\n";
    }

    // Debugging information for use in R.
    if (verbosity > 3) {  // only at very high verbosities
        output += "library(ggplot2)\n";
        output += "library(readr)\n";
        output += "dot.plot.1 <- function(all_points) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size, alpha=p_alpha_1, shape=19) + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
        output += "dot.plot.2 <- function(all_points, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size, alpha=p_alpha_2, shape=19) + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=p_size, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
        output += "dot.plot.3 <- function(all_points, filtered_data, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=p_size, alpha=p_alpha_2, shape=19) + geom_point(data=filtered_data,  aes(x=X1, y=X2), size=p_size, alpha=1, shape=19, colour=\"green\") + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=p_size, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0), limits = x_limits) + scale_y_continuous(expand = c(0, 0), limits = y_limits)}\n";
        output += "p_size <- 0.1\n";
        output += "p_alpha_1 <- 0.1\n";
        output += "p_alpha_2 <- 0.02\n";
    }

    // For each minimap alignment we find the appropriate part of the reference sequence.
    RefRangeMap refRanges;
    for (size_t i = 0; i < minimapAlignments.size(); ++i) {
        std::string minimapStr = minimapAlignments[i];
        std::vector<std::string> minimapStrParts = splitString(minimapStr, ',');

//        std::cout << "`" << minimapStr << "`" << std::endl << std::flush;  // TEMP
        int readStart = std::stoi(minimapStrParts[0]);
        int readEnd = std::stoi(minimapStrParts[1]);
        char readStrand = minimapStrParts[2][0];
        bool posStrand = readStrand == '+';

        std::string refName = minimapStrParts[3];
        int refStart = std::stoi(minimapStrParts[4]);
        int refEnd = std::stoi(minimapStrParts[5]);
        std::string & refSeq = refSeqs->at(refName);

//        std::cout << readStart << "  " << readEnd << "  " << readStrand << "  ";  // TEMP
//        std::cout << refName << "  " << refStart << "  " << refEnd << std::endl;  // TEMP

        StartEndRange refRange = getRefRange(refStart, refEnd, refSeq.length(),
                                             readStart, readEnd, readLength, posStrand);
//        std::cout << refRange.first << "  " << refRange.second << std::endl;  // TEMP

        // The first time we see the reference/strand, initialise it with an empty vector.
        std::string refNameAndStrand = refName + readStrand;
        if (refRanges.find(refNameAndStrand) == refRanges.end())
            refRanges[refNameAndStrand] = std::vector<StartEndRange>();

        refRanges[refNameAndStrand].push_back(refRange);
    }

    // Simplify the reference ranges by combining overlapping ranges.
    RefRangeMap simplifiedRefRanges;
    for(auto const& r : refRanges) {
        std::string refNameAndStrand = r.first;
//        std::cout << refNameAndStrand << std::endl;  // TEMP
        std::vector<StartEndRange> ranges = r.second;
//        for (int i = 0; i < ranges.size(); ++i)  // TEMP
//            std::cout << "(" << ranges[i].first << ","<< ranges[i].second << ") ";  // TEMP
//        std::cout << std::endl;  // TEMP
        std::vector<StartEndRange> simplifiedRanges = simplifyRanges(ranges);
        simplifiedRefRanges[refNameAndStrand] = simplifiedRanges;
//        for (int i = 0; i < simplifiedRanges.size(); ++i)  // TEMP
//            std::cout << "(" << simplifiedRanges[i].first << ","<< simplifiedRanges[i].second << ") ";  // TEMP
//        std::cout << std::endl;  // TEMP
    }

    // Make a new KmerPositions object for the read. We'll actually add positions later as
    // necessary (because we may not need both the positive strand or the negative strand).
    KmerPositions readKmerPositions;
    bool posPositions = false, negPositions = false;

    // Align to each reference range.
    for(auto r : simplifiedRefRanges) {
        std::string refName = r.first;
        char readStrand = refName.back();
        bool posStrand = readStrand == '+';
        refName.pop_back();
        std::string & refSeq = refSeqs->at(refName);
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
        for (auto range : ranges) {
            std::vector<ScoredAlignment *> a =
                alignReadToReferenceRange(refSeqs, refName, range, refSeq.length(), readName,
                                          readStrand, kmerPositions, kSize, readSeq, matchScore,
                                          mismatchScore, gapOpenScore, gapExtensionScore,
                                          sensitivityLevel, verbosity, output);
            returnedAlignments.insert(returnedAlignments.end(), a.begin(), a.end());
        }
    }

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    for (auto alignment : returnedAlignments) {
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
    std::vector<ScoredAlignment *> alignments;

    int refStart = refRange.first;
    int refEnd = refRange.second;
    int readLen = readSeq->length();

//    std::cout << "(" << refName << "," << readStrand << "," << refStart << "," << refEnd << ") ";  // TEMP
//    std::cout << std::flush; // TEMP
    std::string trimmedRefSeq = refSeqs->at(refName).substr(refStart, refEnd-refStart);
    int trimmedRefLen = trimmedRefSeq.length();

    // Find all common k-mer positions.
    std::vector<CommonKmer> commonKmers;
    int maxI = trimmedRefLen - kSize + 1;
    for (int i = 0; i < maxI; ++i) {
        std::string refKmer = trimmedRefSeq.substr(i, kSize);
        if (kmerPositions->find(refKmer) != kmerPositions->end() ) {  // if k-mer is in the read
            std::vector<int> & readPositions = kmerPositions->at(refKmer);
            for (size_t j = 0; j < readPositions.size(); ++j)
                commonKmers.emplace_back(readPositions[j], i);
        }
    }

    // Debugging information for use in R.
    if (verbosity > 3) {  // only at very high verbosities
        std::ofstream allPointsFile;
        std::string filename = readName + readStrand + "_" + refName + "_all_points.tsv";
        allPointsFile.open(filename);
        for (auto k : commonKmers)
            allPointsFile << k.m_hPosition << "\t" << k.m_vPosition << "\n";
        allPointsFile.close();
        output += "all.points <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
    }

    PointCloud cloud;
    addKmerPointsToNanoflann(cloud, commonKmers);
    my_kd_tree_t index(2, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    index.buildIndex();

    double highestDensityScore = 0.0;
    Point highestDensityPoint = getHighestDensityPoint(100, cloud, index, trimmedRefSeq, readSeq,
                                                       &highestDensityScore);

    Point p = highestDensityPoint;
    std::vector<Point> traceDots;
    traceDots.push_back(p);

    int smallLineTracingStepSize = 250;
    int smallSearchRadius = 500;

    // If the line is 'lost' then we will switch to larger steps to increase our change of
    //'finding' it again.
    int largeLineTracingStepSize = 500;
    int largeSearchRadius = 1000;

    // Start the point collection using points around the starting point.
    std::unordered_set<Point> pointSet;
    std::vector<Point> nearbyPoints = radiusSearchAroundPoint(p, smallSearchRadius, cloud, index);
    for (auto nearbyPoint : nearbyPoints) {
        if (closeToDiagonal(p, nearbyPoint))
            pointSet.insert(nearbyPoint);
    }

    int smallestTraceLineX = p.x, largestTraceLineX = p.x;
    int smallestTraceLineY = p.y, largestTraceLineY = p.y;

    // Trace the line forward then backward.
    int directions[2] = {1, -1};
    for (auto direction : directions) {
        p = highestDensityPoint;
        int lineTracingStepSize = smallLineTracingStepSize;
        int searchRadius = smallSearchRadius;
        int maxX = readLen;
        int maxY = trimmedRefSeq.length();
        bool failed;
        while (true) {
            Point newP = p;
            int step = direction * lineTracingStepSize;
            newP.x += step;
            newP.y += step;
            if (newP.x > maxX || newP.y > maxY)
                p = newP;
            else if (newP.x < 0 || newP.y < 0)
                p = newP;
            else {
                p = getHighestDensityPointNearPoint(lineTracingStepSize, newP, cloud, index,
                                                    highestDensityScore, &failed);
                if (failed) {
                    lineTracingStepSize = largeLineTracingStepSize;
                    searchRadius = largeSearchRadius;
                }
                else {
                    lineTracingStepSize = smallLineTracingStepSize;
                    searchRadius = smallSearchRadius;
                }
            }
            if (p.x == -1 || p.y == -1)
                p = newP;

            traceDots.push_back(p);

            // Points near the trace point get added to the point set.
            nearbyPoints = radiusSearchAroundPoint(p, searchRadius, cloud, index);
            for (auto nearbyPoint : nearbyPoints) {
                if (closeToDiagonal(p, nearbyPoint))
                    pointSet.insert(nearbyPoint);
            }

            smallestTraceLineX = std::min(p.x, smallestTraceLineX);
            smallestTraceLineY = std::min(p.y, smallestTraceLineY);
            largestTraceLineX = std::max(p.x, largestTraceLineX);
            largestTraceLineY = std::max(p.y, largestTraceLineY);

            if (direction == 1 && (p.x > maxX || p.y > maxY))
                break;
            if (direction == -1 && (p.x < 0 || p.y < 0))
                break;
        }
    }

    if (verbosity > 3) {  // only at very high verbosities
        std::ofstream traceDotsFile;
        std::string filename = readName + readStrand + "_" + refName + "_trace_dots.tsv";
        traceDotsFile.open(filename);
        for (auto d : traceDots)
            traceDotsFile << d.x << "\t" << d.y << "\n";
        traceDotsFile.close();
        output += "trace.dots <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
        
        std::ofstream filteredDataFile;
        filename = readName + readStrand + "_" + refName + "_filtered_data.tsv";
        filteredDataFile.open(filename);
        for (auto d : pointSet)
            filteredDataFile << d.x << "\t" << d.y << "\n";
        filteredDataFile.close();
        output += "filtered.data <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
    }

    // Now we can do a Seqan alignment around the points we've collected!

    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed> TSeedSet;

    String<TSeed> seeds;
    for (auto p : pointSet)
        appendValue(seeds, TSeed(p.x, p.y, kSize));

    TSeedSet seedSet;
    for (unsigned i = 0; i < length(seeds); ++i) {
        if (!addSeed(seedSet, seeds[i], 2, Merge()))
            addSeed(seedSet, seeds[i], Single());
    }

    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());

    if (verbosity > 3) {  // only at very high verbosities
        std::ofstream chainedSeedsFile;
        std::string filename = readName + readStrand + "_" + refName + "_chained_seeds.tsv";
        chainedSeedsFile.open(filename);
        for (unsigned i = 0; i < length(seedChain); ++i)
            chainedSeedsFile << beginPositionH(seedChain[i]) << "\t" << beginPositionV(seedChain[i]) << "\n";
        chainedSeedsFile.close();
        output += "chained.seeds <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
        output += "x_limits <- c(0, max(all.points[,1]))\n";
        output += "y_limits <- c(0, max(all.points[,2]))\n";
        output += "dot.plot.1(all.points)\n";
        output += "dot.plot.2(all.points, trace.dots)\n";
        output += "dot.plot.3(all.points, filtered.data, trace.dots)\n";
        output += "dot.plot.3(all.points, chained.seeds, trace.dots)\n";
    }

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), *readSeq);
    assignSource(row(alignment, 1), trimmedRefSeq);
    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    int bandSize = LEVEL_0_BAND_SIZE;
    if (sensitivityLevel == 1)
        bandSize = LEVEL_1_BAND_SIZE;
    else if (sensitivityLevel == 2)
        bandSize = LEVEL_2_BAND_SIZE;
    else if (sensitivityLevel == 3)
        bandSize = LEVEL_3_BAND_SIZE;

    // First try a fast alignment with a small band.
    ScoredAlignment * sgAlignment;
    try {
        bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);
        std::string signedReadName = readName + readStrand;
        sgAlignment = new ScoredAlignment(alignment, signedReadName, refName, readLen, refLen,
                                          refStart, startTime, bandSize, false, false, false,
                                          scoringScheme);
    }
    catch (...) {
        sgAlignment = 0;
    }

    alignments.push_back(sgAlignment);
    return alignments;
}


void addKmerPointsToNanoflann(PointCloud & cloud, std::vector<CommonKmer> & commonKmers) {
    cloud.pts.resize(commonKmers.size());
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        cloud.pts[i].x = commonKmers[i].m_hPosition;
        cloud.pts[i].y = commonKmers[i].m_vPosition;
    }
}


std::vector<Point> radiusSearchAroundPoint(Point point, int radius, PointCloud & cloud,
                                           my_kd_tree_t & index) {
    std::vector<Point> points;
    nanoflann::SearchParams params;
    std::vector<std::pair<size_t,int> > ret_matches;
    const int query_pt[2] = {point.x, point.y};
    index.radiusSearch(query_pt, radius, ret_matches, params);
    for (auto i : ret_matches)
        points.push_back(cloud.pts[i.first]);
    return points;
}

std::vector<Point> getPointsInHighestDensityRegion(int searchRadius, std::string & trimmedRefSeq,
                                                   std::string * readSeq, PointCloud & cloud,
                                                   my_kd_tree_t & index) {

    int xStepCount = int(ceil(readSeq->length() / double(searchRadius)));
    int yStepCount = int(ceil(trimmedRefSeq.length() / double(searchRadius)));
    double xStepSize = double(readSeq->length()) / xStepCount;
    double yStepSize = double(trimmedRefSeq.length()) / yStepCount;

//    std::cout << std::endl;
//    std::cout << "xStepCount=" << xStepCount << std::endl;
//    std::cout << "yStepCount=" << yStepCount << std::endl;
//    std::cout << "xStepSize=" << xStepSize << std::endl;
//    std::cout << "yStepSize=" << yStepSize << std::endl << std::flush;

    nanoflann::SearchParams params;
    double highestDensity = 0.0;
    std::vector<Point> pointsInHighestDensity;

    for (int i = 0; i <= xStepCount; ++i) {
        int xCentre = int(0.5 + i * xStepSize);

        for (int j = 0; j <= yStepCount; ++j) {
            int yCentre = int(0.5 + j * yStepSize);

//            std::cout << "xCentre=" << xCentre << ", yCentre=" << yCentre << std::endl << std::flush;

            const int query_pt[2] = {xCentre, yCentre};

            std::vector<std::pair<size_t,int> > ret_matches;
            const size_t nMatches = index.radiusSearch(query_pt, searchRadius, ret_matches, params);
            double density = double(nMatches);

//            std::cout << "unadjusted density=" << density << std::endl << std::flush;

            // If the search region is on the edge, increase the density (because the region has
            // less area). It would be technically correct to double the density on the edges,
            // but this biases the density peak towards edges, so we only use 1.5 instead.
            if (i == 0 || i == xStepCount) density *= 1.5;
            if (j == 0 || j == yStepCount) density *= 1.5;

//            std::cout << "adjusted density=" << density << std::endl << std::flush;

            if (density > highestDensity) {
//                std::cout << "NEW BEST!" << std::endl << std::flush;
                highestDensity = density;
                pointsInHighestDensity.clear();
                for (auto k : ret_matches)
                    pointsInHighestDensity.push_back(cloud.pts[k.first]);
            }
        }
    }
    return pointsInHighestDensity;
}

Point getHighestDensityPoint(int densityRadius, PointCloud & cloud, my_kd_tree_t & index,
                             std::string & trimmedRefSeq, std::string * readSeq,
                             double * highestDensityScore) {

    std::vector<Point> points = getPointsInHighestDensityRegion(densityRadius * 2, trimmedRefSeq,
                                                                readSeq, cloud, index);
    Point highestDensityPoint = points[0];
    *highestDensityScore = 0.0;

    for (auto point : points) {
        double densityScore = getPointDensityScore(densityRadius, point, cloud, index);
        if (densityScore > *highestDensityScore) {
            *highestDensityScore = densityScore;
            highestDensityPoint = point;
        }
    }
    return highestDensityPoint;
}

Point getHighestDensityPointNearPoint(int densityRadius, Point centre, PointCloud & cloud,
                                      my_kd_tree_t & index, double highestDensityScore,
                                      bool * failed) {
    std::vector<Point> points = radiusSearchAroundPoint(centre, densityRadius, cloud, index);
    if (points.size() == 0)
        return {-1, -1};
    Point highestDensityPoint = centre;
    *failed = true;
    double bestDensityScore = highestDensityScore / 10.0;

    for (auto point : points) {
        double densityScore = getPointDensityScore(densityRadius, point, cloud, index);

        // Boost the density score for points near the centre.
        int distanceFromCentre = abs(point.x - centre.x) + abs(point.y - centre.y);
        densityScore *= (1.0 + ((densityRadius - distanceFromCentre) / densityRadius));

        if (densityScore > bestDensityScore) {
            bestDensityScore = densityScore;
            highestDensityPoint = point;
            *failed = false;
        }
    }

//    std::cout << "Starting point: " << centre.x << "," << centre.y << "\n";
//    std::cout << "Highest density point: " << highestDensityPoint.x << "," << highestDensityPoint.y << "\n";
//    std::cout << "Density score: " << bestDensityScore << "\n";
//    std::cout << "\n";

    return highestDensityPoint;
}


double getPointDensityScore(int densityRadius, Point p, PointCloud & cloud, my_kd_tree_t & index) {
    std::vector<Point> neighbourPoints = radiusSearchAroundPoint(p, densityRadius, cloud, index);
    double densityScore = 0.0;
    for (auto neighbourPoint : neighbourPoints) {
        int xDiff = neighbourPoint.x - p.x;
        int yDiff = neighbourPoint.y - p.y;
        if (xDiff + yDiff > 0)
            densityScore += 1.0 / (abs(xDiff-yDiff) + 1.0);
    }
    return densityScore;
}


double fractionOfReadAligned(std::vector<ScoredAlignment *> & alignments) {
    if (alignments.size() == 0)
        return true;
    std::vector<std::pair<int, int> > ranges;
    for (size_t i = 0; i < alignments.size(); ++i) {
        ScoredAlignment * alignment = alignments[i];
        int start, end;
        if (alignment->isRevComp()) {
            start = alignment->m_readLength - alignment->m_readEndPos;
            end = alignment->m_readLength - alignment->m_readStartPos;
        }
        else {
            start = alignment->m_readStartPos;
            end = alignment->m_readEndPos;
        }
        ranges.push_back(std::pair<int, int>(start, end));
    }
    std::vector<std::pair<int, int> > simplifiedRanges = simplifyRanges(ranges);
    int alignedLength = 0;
    for (size_t i = 0; i < simplifiedRanges.size(); ++i)
        alignedLength += simplifiedRanges[i].second - simplifiedRanges[i].first;
    return double(alignedLength) / alignments[0]->m_readLength;
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
    std::sort(ranges.begin(),ranges.end());
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


bool closeToDiagonal(Point p1, Point p2) {
    int xDiff = p1.x - p2.x;
    int yDiff = p1.y - p2.y;
    double slope = 0.0;
    if (xDiff != 0)
        slope = double(yDiff) / double(xDiff);
    if (xDiff == 0 && yDiff == 0)
        slope = 1.0;
    return (slope > 0.6667 && slope < 1.5);
}
