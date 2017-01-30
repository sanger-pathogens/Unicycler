# This makefile will build the C++ components of Unicycler.

# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the *.so file, which is required to run the aligner)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXXFLAGS can be overridden by the user.
CXXFLAGS    ?= -Wall -Wextra -pedantic -march=native

# These flags are required for the build to work.
FLAGS        = -std=c++14 -Iunicycler/include -fPIC
LDFLAGS      = -shared -lz

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = unicycler/cpp_functions.so
SHELL        = /bin/sh
SOURCES      = $(shell find unicycler -name "*.cpp")
HEADERS      = $(shell find unicycler -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)

# Linux needs '-soname' while Mac needs '-install_name'
PLATFORM     = $(shell uname)
ifeq ($(PLATFORM), Darwin)
SONAME       = -install_name
else
SONAME       = -soname
endif

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(LDFLAGS) -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
