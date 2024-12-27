#ifndef MATCH_TOOLS_HPP
#define MATCH_TOOLS_HPP

#include "Definitions.hpp"

#include "TString.h"

int MatchIdx(LorentzVectorF_t const& parton, VecLVF_t const& jets);
TString MakeTrueLabel(VecLVF_t const& gen, VecLVF_t const& reco);
inline bool IsTrue(TString const& true_label, TString const& label) { return true_label == label; }

#endif