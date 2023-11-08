#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"

void save_1d_dist(TH1F* dist, std::string const& name, std::string const& title);
void save_2d_dist(TH2F* dist, std::string const& name, std::string const& title_x, std::string const& title_y);

#endif