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
#define LEVEL_3_MIN_LINE_TRACE_COUNT 4

// The max line trace count controls how many line tracings will possibly be tried for each
// alignment, assuming previous line tracings looked bad. E.g. if set to 10, up to 10 line tracings
// will be tried, but quite possibly fewer will actually be done because the code will stop trying
// when a line looks good.
#define LEVEL_0_MAX_LINE_TRACE_COUNT 8
#define LEVEL_1_MAX_LINE_TRACE_COUNT 16
#define LEVEL_2_MAX_LINE_TRACE_COUNT 32
#define LEVEL_3_MAX_LINE_TRACE_COUNT 64
