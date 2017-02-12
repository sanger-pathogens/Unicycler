// Copyright 2017 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Unicycler

// This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Unicycler is
// distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Unicycler. If not, see <http://www.gnu.org/licenses/>.

#define LEVEL_0_MINIMAP_KMER_SIZE 15
#define LEVEL_1_MINIMAP_KMER_SIZE 15
#define LEVEL_2_MINIMAP_KMER_SIZE 14
#define LEVEL_3_MINIMAP_KMER_SIZE 13

#define LEVEL_0_KMER_SIZE 10
#define LEVEL_1_KMER_SIZE 10
#define LEVEL_2_KMER_SIZE 9
#define LEVEL_3_KMER_SIZE 8

#define LEVEL_0_BAND_SIZE 25
#define LEVEL_1_BAND_SIZE 50
#define LEVEL_2_BAND_SIZE 75
#define LEVEL_3_BAND_SIZE 100

// The min line trace count controls how many line tracings will definitely be tried for each
// alignment. If more than 1, then additional line tracings will be tried, even if the first one
// looked good.
#define LEVEL_0_MIN_LINE_TRACE_COUNT 1
#define LEVEL_1_MIN_LINE_TRACE_COUNT 1
#define LEVEL_2_MIN_LINE_TRACE_COUNT 2
#define LEVEL_3_MIN_LINE_TRACE_COUNT 3

// The max line trace count controls how many line tracings will possibly be tried for each
// alignment, assuming previous line tracings looked bad. E.g. if set to 10, up to 10 line tracings
// will be tried, but quite possibly fewer will actually be done because the code will stop trying
// when a line looks good.
#define LEVEL_0_MAX_LINE_TRACE_COUNT 4
#define LEVEL_1_MAX_LINE_TRACE_COUNT 8
#define LEVEL_2_MAX_LINE_TRACE_COUNT 12
#define LEVEL_3_MAX_LINE_TRACE_COUNT 16

// When beginning line tracing, the location with the highest density is found by checking many
// regions of the space. This parameter controls the size of those regions.
#define LINE_TRACING_START_POINT_SEARCH_RADIUS 100

// Common k-mer points this close to the trace line will be collected into the point set that is
// given to Seqan for global chaining and then banded alignment.
#define TRACE_LINE_COLLECTION_DISTANCE 20.0

// The trace line points will be this far apart.
#define TRACE_LINE_STEP_DISTANCE 500

// This is how much trace line points change when we're searching for the best position.
#define TRACE_LINE_MUTATION_SIZE 5

// This is the score a line segment will get it if it exactly follows points for its whole length.
// Increase this value to make the line finding more strongly prefer point proximity (as opposed
// to slope).
#define MAX_POINTS_SCORE 10000.0

// This is the maximum score a line segment can get for its slope (achieved with a slope of 1).
// Increase this value to make the line finding more strongly prefer slope (as opposed to point
// proximity).
#define MAX_SLOPE_SCORE 1.0

// If a line segment's slope (or its reciprocal) is less than this, it gets a slope score of 0.
#define MIN_ACCEPTABLE_LINE_SEGMENT_SLOPE 0.5

// If a Seqan seed chain contains a gap with an area larger than this, then we don't go ahead with
// the alignment (because it would take too long and probably not be good anyway).
#define MAX_BANDED_ALIGNMENT_GAP_AREA 100000000
