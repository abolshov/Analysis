CXX=clang++
CXXFLAGS= -c -O2 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

analysis.o: analysis.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Analyzer.o: Analyzer.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Storage.o: Storage.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Estimator.o: Estimator.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimatorUtils.o: EstimatorUtils.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimatorTools.o: EstimatorTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

MatchingTools.o: MatchingTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

SelectionUtils.o: SelectionUtils.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

HistManager.o: HistManager.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimationRecorder.o: EstimationRecorder.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

analysis: analysis.o Analyzer.o Storage.o Estimator.o EstimatorUtils.o EstimatorTools.o MatchingTools.o HistManager.o SelectionUtils.o EstimationRecorder.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm analysis
	rm *.o