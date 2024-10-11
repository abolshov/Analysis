#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include <memory>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "Constants.hpp"

using UHist1D = std::unique_ptr<TH1F>;
using UHist2D = std::unique_ptr<TH2F>;

struct EstimatorInput
{
    EstimatorInput(std::vector<TLorentzVector>&& p4, std::string const& pdf_file);

    std::vector<TLorentzVector> p4;
    std::vector<UHist1D> pdf1d;
    std::vector<UHist2D> pdf2d;
};

#endif