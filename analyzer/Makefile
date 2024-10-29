CXX=clang++
CXXFLAGS= -c -O2 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

analysis.o: analysis.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Analyzer.o: Analyzer.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Event.o: Event.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Objects.o: Objects.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Validator.o: Validator.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimatorTools.o: EstimatorTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

HistManager.o: HistManager.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

analysis: analysis.o Analyzer.o Event.o Objects.o Validator.o EstimatorTools.o HistManager.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm analysis
	rm *.o