#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include "TRandom3.h"
#include "Definitions.hpp"

LorentzVectorF_t SamplePNetResCorr(LorentzVectorF_t const& jet, std::unique_ptr<TRandom3>& prg, Float_t resolution);
std::pair<Float_t, Float_t> ComputeJetResc(LorentzVectorF_t const& b1, LorentzVectorF_t const& b2, UHist_t<TH1F>& pdf, Float_t mass);
std::optional<LorentzVectorF_t> NuFromOnshellW(Float_t eta, Float_t phi, Float_t mw, LorentzVectorF_t const& lep_onshell);
std::optional<LorentzVectorF_t> NuFromOffshellW(LorentzVectorF_t const& lep1, 
                                                LorentzVectorF_t const& lep2, 
                                                LorentzVectorF_t const& nu1,
                                                LorentzVectorF_t const& met,
                                                int control, 
                                                Float_t mh);

#endif