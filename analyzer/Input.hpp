#ifndef INPUT_HPP
#define INPUT_HPP

#include <vector>
#include <memory>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "Constants.hpp"
#include "Event.hpp"

using UHist1D = std::unique_ptr<TH1F>;
using UHist2D = std::unique_ptr<TH2F>;

using Vec1D_t = std::vector<UHist1D>;
using Vec2D_t = std::vector<UHist2D>;

struct EstimatorInput
{
    EstimatorInput() = default;
    EstimatorInput(std::vector<TLorentzVector>&& p4, std::string const& pdf_file);
    EstimatorInput(std::vector<TLorentzVector>&& p4, std::vector<UHist1D>&& vec_pdf_1d, std::vector<UHist2D>&& vec_pdf_2d);

    std::vector<TLorentzVector> p4; 
    std::vector<UHist1D> pdf1d;
    std::vector<UHist2D> pdf2d;
};

struct ValidatorInput
{
    ValidatorInput(Event const& event, Vec1D_t&& pdf_1d, Vec2D_t&& pdf_2d);
    Vec1D_t pdf1d;
    Vec2D_t pdf2d;
    std::vector<TLorentzVector> gen_truth_p4; // p4 of gen level quarks, lepton and met
    std::vector<TLorentzVector> reco_jet_p4; // p4 of all selected reco jets
    std::vector<TLorentzVector> nu;
    std::vector<double> jet_PNet_resolutions; // pt resolutions provided by pnet
    std::vector<double> jet_PNet_corrections; // central rescaling value provided by pnet
    TLorentzVector recoMET;
};

#endif