#ifndef MATCH_TOOLS_HPP
#define MATCH_TOOLS_HPP

#include "Definitions.hpp"
#include "JetCombination.hpp"
#include "Storage.hpp"

#include "TString.h"

int MatchIdx(LorentzVectorF_t const& parton, VecLVF_t const& jets);

// to ensure compatibility with how label of any combination is formed
// will organize indices of matched reco b jets lexicographically
// will organize indices of matched reco light jets lexicographically
TString MakeTrueLabel(VecLVF_t const& gen, VecLVF_t const& reco);
inline bool IsTrue(TString const& true_label, TString const& label) { return true_label == label; }

JetComb FindMatch(VecLVF_t const& quarks, VecLVF_t const& jets, Channel ch);

#endif