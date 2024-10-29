#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP

#include "Event.hpp"
#include "Definitions.hpp"
#include "TFile.h"

class Validator 
{
    public:
    void FillVariables(Event const& event);
    inline void Set1dPDF(HistVec1d_t&& pdf1d) { pdf_1d = std::move(pdf1d); }
    inline void Set2dPDF(HistVec2d_t&& pdf2d) { pdf_2d = std::move(pdf2d); }
    inline void SetTruthIndex(std::map<std::string, size_t> const& index) { truth_index = index; }
    inline void SetPDF1dIndex(std::map<std::string, size_t> const& index) { pdf1d_index = index; }
    inline void SetPDF2dIndex(std::map<std::string, size_t> const& index) { pdf2d_index = index; }

    void ValidateBJetCorr();
    void Print();

    private:
    std::map<std::string, size_t> truth_index;
    std::map<std::string, size_t> pdf1d_index;
    std::map<std::string, size_t> pdf2d_index;

    HistVec1d_t pdf_1d;
    HistVec2d_t pdf_2d;

    VecLVF_t gen_truth_p4; // p4 of gen level quarks, lepton and met
    VecLVF_t reco_jet_p4; // p4 of all selected reco jets
    VecLVF_t nu;

    std::vector<Float_t> jet_PNet_resolutions; // pt resolutions provided by pnet
    std::vector<Float_t> jet_PNet_corrections; // central rescaling value provided by pnet
    
    LorentzVectorF_t recoMET; 
};

#endif