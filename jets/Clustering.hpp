#ifndef CLUST_HPP
#define CLUST_HPP

#include <vector>
#include <algorithm>
#include <set>

#include "TROOT.h"

#include "Index.hpp"

static const std::set<Int_t> incoming{ 12, 13, 21, 31, 41, 41, 45, 46, 53, 61, 158 };

std::vector<Int_t> GetFinalParticles(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart); // returns positions of final state particles
std::vector<Bool_t> PossibleJetConstituents(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart, GenPartIndex const& idx); // marks all mothers as NOT final state (bc they decay)

inline Bool_t IsIncoming(Int_t status) { return incoming.find(status) != incoming.end(); }

#endif