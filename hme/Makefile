CXX=clang++
CXXFLAGS= -c -O3 -std=c++17 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

hme.o: hme.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimatorTools.o: EstimatorTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

HistManager.o: HistManager.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

MatchingTools.o: MatchingTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

hme: hme.o EstimatorTools.o HistManager.o MatchingTools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm hme
	rm *.o
