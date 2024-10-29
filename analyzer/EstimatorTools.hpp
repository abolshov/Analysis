#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include "TRandom3.h"

#include "Definitions.hpp"
#include "Constants.hpp"

// computes rescaling factors for jets using rescale pdf and mass of particle->jj
// has side effect: may swap j1 and j2
OptionalPair JetRescFact(LorentzVectorF_t& j1, LorentzVectorF_t& j2, UHist1d_t& pdf, double mass);

#endif