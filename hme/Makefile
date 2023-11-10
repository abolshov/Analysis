ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS)
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS)

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .C
.PREFIXES: ./obj/

#----------------------------------------------------#

hme: hme.cpp tools.cpp Utils.cpp
	g++ -O2 -Wall -Wextra -Wpedantic $^ -o $@ `root-config --cflags --glibs ` -lSpectrum

compute: compute_pdfs.cpp Utils.cpp Plotting.cpp
	g++ -O2 -Wall -Wextra -Wpedantic $^ -o $@ `root-config --cflags --glibs ` -lSpectrum 

test: test.cpp
	g++ -O3 -Wall -Wextra -Wpedantic $^ -o $@ `root-config --cflags --glibs ` -lSpectrum 
	
clean:
	rm -f hme
	rm -f *~
	rm *.o 
