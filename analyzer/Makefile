CXX=clang++
CXXFLAGS= -c -O2 -Wall -Wextra -pedantic `root-config --cflags `
LDFLAGS= `root-config --glibs ` -lSpectrum

analysis.o: analysis.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Analyzer.o: Analyzer.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Storage.o: Storage.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

analysis: analysis.o Analyzer.o Storage.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean: 
	rm analysis
	rm *.o