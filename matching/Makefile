CXX=g++
CXXFLAGS= -c -O3 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

matching.o: matching.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

MatchingTools.o: MatchingTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

matching: matching.o MatchingTools.o
	g++ $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm matching
	rm *.o
