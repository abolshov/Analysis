CXX=clang++
CXXFLAGS= -c -O3 -Wall -Wextra -Wpedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

matching.o: matching.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

MatchingTools.o: MatchingTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

HistManager.o: HistManager.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

matching: matching.o MatchingTools.o HistManager.o
	g++ $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm matching
	rm *.o
