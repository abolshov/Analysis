#ifndef CLUST_HPP
#define CLUST_HPP

#include <vector>
#include <algorithm>
#include <set>

#include "TROOT.h"
#include "TH1.h"

#include "Index.hpp"

static const std::set<Int_t> incoming{ 12, 13, 21, 31, 41, 41, 45, 46, 53, 61, 158 };

std::vector<Int_t> GetFinalParticles(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart); // returns positions of final state particles
std::vector<Bool_t> PossibleJetConstituents(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart, GenPartIndex const& idx); // marks all mothers as NOT final state (bc they decay)
// std::vector<Int_t> GetDaughtersOf(Int_t mother_pdgId, Int_t mother_idx);

inline Bool_t IsIncoming(Int_t status) { return incoming.find(status) != incoming.end(); }
Bool_t IsFinal(Int_t pos, Int_t const* motherIdxs, Int_t nGenPart);
Bool_t DecayProductOf(Int_t this_part, Int_t pos, Int_t const* motherIdxs); // is this_part product of decay of particle at pos

void AnalyzeJets(std::unique_ptr<TH1I> const& unused_cand, std::unique_ptr<TH1I> const& bad_jets, PtEtaPhiMArray const& genPart, PtEtaPhiMArray const& genJet, std::vector<Bool_t> candidates, Int_t n_unused);

std::vector<Overlap> FindOverlaps(PtEtaPhiMArray const& gen_jets, PtEtaPhiMArray const& gen_parts, std::vector<Bool_t> const& finals);
std::pair<Int_t, Int_t> FindBestMatch(Overlap const& ov, PtEtaPhiMArray const& gen_jets, PtEtaPhiMArray const& gen_parts);

#endif