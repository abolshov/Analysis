CXX=clang++
CXXFLAGS= -c -O2 -std=c++17 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

reco_hme.o: reco_hme.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

EstimatorTools.o: EstimatorTools.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

HistManager.o: HistManager.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

reco_hme: reco_hme.o EstimatorTools.o HistManager.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm hme
	rm *.o

compute_pdf: compute_pdf.cpp
	g++ -Wall -Wextra -O2 -o compute_pdf compute_pdf.cpp `root-config --cflags --glibs ` -lSpectrum