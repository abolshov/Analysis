compute_pdf: compute_pdf.cpp
	g++ -Wall -Wextra -O2 -o compute_pdf compute_pdf.cpp `root-config --cflags --glibs ` -lSpectrum

reco_hme: reco_hme.cpp HistManager.cpp	
	g++ -Wall -Wextra -O2 -o reco_hme reco_hme.cpp HistManager.cpp `root-config --cflags --glibs ` -lSpectrum