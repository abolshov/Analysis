CXX=g++
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

analysis: analysis.o Analyzer.o Storage.o Estimator.o EstimatorUtils.o EstimatorTools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm analysis
	rm *.o