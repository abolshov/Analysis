CXX=g++
CXXFLAGS= -c -O3 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

gather.o: gather.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

MatchingTools.o: MatchingTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

gather: gather.o MatchingTools.o
	g++ $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm gather
	rm *.o
