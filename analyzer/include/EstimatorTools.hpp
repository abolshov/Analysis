#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include "TRandom3.h"
#include "Definitions.hpp"

LorentzVectorF_t SamplePNetResCorr(LorentzVectorF_t const& jet, std::unique_ptr<TRandom3>& prg, Float_t resolution);
OptPairF_t ComputeJetResc(LorentzVectorF_t const& p1, LorentzVectorF_t const& p2, UHist_t<TH1F>& pdf, Float_t mass);
OptLVecF_t NuFromOnshellW(Float_t eta, Float_t phi, Float_t mw, LorentzVectorF_t const& lep_onshell);
OptLVecF_t NuFromOffshellW(LorentzVectorF_t const& lep1, LorentzVectorF_t const& lep2, LorentzVectorF_t const& nu1, LorentzVectorF_t const& met, int control, Float_t mh);
OptLVecF_t NuFromH(LorentzVectorF_t const& jet1, LorentzVectorF_t const& jet2, LorentzVectorF_t const& lep, LorentzVectorF_t const& met, bool add_deta, Float_t mh);
OptLVecF_t NuFromW(LorentzVectorF_t const& lep, LorentzVectorF_t const& met, bool add_deta, Float_t mw);

inline Float_t mT(LorentzVectorF_t p)
{
    return std::sqrt(p.M2() + p.Pt()*p.Pt());
}

OptLVecFPair_t NuFromConstraints(LorentzVectorF_t const& jet1, LorentzVectorF_t const& jet2, 
                                LorentzVectorF_t const& lep, LorentzVectorF_t const& met, 
                                Float_t mw, Float_t mh);

#endif