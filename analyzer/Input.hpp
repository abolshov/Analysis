#ifndef INPUT_HPP
#define INPUT_HPP

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "Constants.hpp"
#include "Event.hpp"
#include "Definitions.hpp"


struct EstimatorInput
{
    EstimatorInput() = default;
    EstimatorInput(Event const& event, HistVec1d_t&& vec_pdf_1d, HistVec2d_t&& vec_pdf_2d);
 
    HistVec1d_t pdf1d;
    HistVec2d_t pdf2d;

    VecLVF_t p4;
};

struct ValidatorInput
{
    ValidatorInput() = default;
    ValidatorInput(Event const& event, HistVec1d_t&& pdf_1d, HistVec2d_t&& pdf_2d);

    HistVec1d_t pdf1d;
    HistVec2d_t pdf2d;

    VecLVF_t gen_truth_p4; // p4 of gen level quarks, lepton and met
    VecLVF_t reco_jet_p4; // p4 of all selected reco jets
    VecLVF_t nu;

    std::vector<Float_t> jet_PNet_resolutions; // pt resolutions provided by pnet
    std::vector<Float_t> jet_PNet_corrections; // central rescaling value provided by pnet
    
    LorentzVectorF_t recoMET;
};

#endif