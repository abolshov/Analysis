#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include "TRandom3.h"
#include "Definitions.hpp"

LorentzVectorF_t SamplePNetResCorr(LorentzVectorF_t const& jet, std::unique_ptr<TRandom3>& prg, Float_t resolution);
std::pair<Float_t, Float_t> ComputeJetResc(LorentzVectorF_t const& p1, LorentzVectorF_t const& p2, UHist_t<TH1F>& pdf, Float_t mass);
std::optional<LorentzVectorF_t> NuFromOnshellW(Float_t eta, Float_t phi, Float_t mw, LorentzVectorF_t const& lep_onshell);
std::optional<LorentzVectorF_t> NuFromOffshellW(LorentzVectorF_t const& lep1, LorentzVectorF_t const& lep2, LorentzVectorF_t const& nu1, LorentzVectorF_t const& met, int control, Float_t mh);
std::optional<LorentzVectorF_t> NuFromHiggsConstr(LorentzVectorF_t const& jet1, LorentzVectorF_t const& jet2, LorentzVectorF_t const& lep, LorentzVectorF_t const& met, int control, Float_t mh);
std::optional<LorentzVectorF_t> NuFromWConstr(LorentzVectorF_t const& lep, LorentzVectorF_t const& met, int control, Float_t mw);

inline Float_t mT(LorentzVectorF_t p)
{
    return std::sqrt(p.M2() + p.Pt()*p.Pt());
}

#endif