CXX=g++
CXXFLAGS= -c -O3 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

analysis.o: analysis.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Analyzer.o: Analyzer.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EventData.o: EventData.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Selector.o: Selector.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

analysis: analysis.o Analyzer.o EventData.o Selector.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm analysis
	rm *.o