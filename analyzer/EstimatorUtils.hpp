#ifndef ESTIMATOR_UTILS_HPP
#define ESTIMATOR_UTILS_HPP

#include <vector>

#include "Definitions.hpp"
#include "Storage.hpp"

#include "TFile.h"
#include "TH1.h"

template <typename T>
UHist_t<T> ReadHist(TFile* file, TString hist_name)
{
    UHist_t<T> res(static_cast<T*>(file->Get<T>(hist_name)));
    res->SetDirectory(nullptr);
    return res;
}

void Get1dPDFs(TFile* fptr, HistVec_t<TH1F>& pdfs);
void Get2dPDFs(TFile* fptr, HistVec_t<TH2F>& pdfs);

inline LorentzVectorF_t GetRecoMET(Storage const& s) { return LorentzVectorF_t(s.reco_met_pt, 0.0, s.reco_met_phi, 0.0); }
inline LorentzVectorF_t GetGenMET(Storage const& s) { return LorentzVectorF_t(s.gen_met_pt, 0.0, s.gen_met_phi, 0.0); }

VecLVF_t GetRecoJetP4(Storage const& s);
VecLVF_t GetRecoLepP4(Storage const& s, Channel ch);

// IMPLEMENT LATER
// VecLVF_t GetGenJetP4(Storage const& s);
// VecLVF_t GetGenLepP4(Storage const& s);

std::vector<Float_t> GetPNetRes(Storage const& s);

#endif