#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include "TRandom3.h"
#include "Definitions.hpp"

LorentzVectorF_t SamplePNetResCorr(LorentzVectorF_t const& jet, std::unique_ptr<TRandom3>& prg, Float_t resolution);
std::pair<Float_t, Float_t> ComputeJetResc(LorentzVectorF_t const& b1, LorentzVectorF_t const& b2, UHist_t<TH1F>& pdf, Float_t mass);

#endif