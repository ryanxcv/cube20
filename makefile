CXXFLAGS = -DHALF -std=c++11 -O3 -Wall -DLEVELCOUNTS -DTHREADS
LIBS = -lpthread

default: all

BINARIES = scramble cubeutil

all: $(BINARIES)

test: $(TESTBINARIES)
	./cubepos_test && ./kocsymm_test && ./phase2prune_test && ./phase1prune_test

scramble: scramble.cpp table.h twophase.h phase1prune.cpp phase1prune.h phase2prune.cpp phase2prune.h kocsymm.cpp kocsymm.h cubepos.cpp cubepos.h
	$(CXX) $(CXXFLAGS) -o $@ scramble.cpp phase1prune.cpp phase2prune.cpp kocsymm.cpp cubepos.cpp $(LIBS)

cubeutil: cubeutil.cpp table.h phase1prune.cpp phase1prune.h kocsymm.cpp kocsymm.h cubepos.cpp cubepos.h
	$(CXX) $(CXXFLAGS) -o $@ cubeutil.cpp phase1prune.cpp kocsymm.cpp cubepos.cpp $(LIBS)
