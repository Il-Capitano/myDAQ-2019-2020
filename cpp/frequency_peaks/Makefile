# compiler parameters
CXX := g++
CXX_FLAGS_DEBUG := -c -g -Og -Wall -Wextra -std=c++17
CXX_FLAGS_RELEASE := -c -O3 -Wall -Wextra -std=c++17
CXX_INCLUDES :=

CC := gcc
CC_FLAGS_DEBUG := -c -Og -g -Wall -Wextra -std=c11
CC_FLAGS_RELEASE := -c -O3 -Wall -Wextra -std=c11
CC_INCLUDES :=

# linker parameters
LD := g++
LD_FLAGS_DEBUG := -g -shared
LD_FLAGS_RELEASE := -shared
LD_LIBDIRS :=
LD_LIBS :=

# debug executable path
DLL_DEBUG := bin/debug/frequency_peaks.dll

# release executable path
DLL_RELEASE := bin/release/frequency_peaks.dll

# all objs that are in BIN_INT
# used for clean-up
ALL_OBJ := $(wildcard bin/debug/int/*.o) $(wildcard bin/release/int/*.o)
ALL_EXE := $(wildcard bin/debug/*.dll) $(wildcard bin/release/*.dll)

# object files and executable are rebuildable
# used for clean-up
REBUILDABLES := $(ALL_EXE) $(ALL_OBJ)

# delete and copy commands
RM := rm
CP := cp

# builds the debug executable
all: debug

# debug build
debug: $(DLL_DEBUG)

# release build
release: $(DLL_RELEASE)

# cleans the binary files
clean:
	$(RM) $(REBUILDABLES)

$(DLL_DEBUG): bin/debug/int/find_peaks.cpp.o
	$(LD) $(LD_FLAGS_DEBUG) $(addprefix -L,$(LD_LIBDIRS)) $(addprefix -l,$(LD_LIBS)) $^ -o $@
	$(CP) $@ ../../labview/spectrum_analysis

$(DLL_RELEASE): bin/release/int/find_peaks.cpp.o
	$(LD) $(LD_FLAGS_RELEASE) $(addprefix -L,$(LD_LIBDIRS)) $(addprefix -l,$(LD_LIBS)) $^ -o $@
	$(CP) $@ ../../labview/spectrum_analysis

bin/debug/int/find_peaks.cpp.o: ./src/find_peaks.cpp
	$(CXX) $(CXX_FLAGS_DEBUG) $(addprefix -I,$(CXX_INCLUDES)) $< -o $@

bin/release/int/find_peaks.cpp.o: ./src/find_peaks.cpp
	$(CXX) $(CXX_FLAGS_RELEASE) $(addprefix -I,$(CXX_INCLUDES)) $< -o $@

