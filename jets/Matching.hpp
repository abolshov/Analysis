#ifndef MATCHING_HPP
#define MATCHING_HPP

#include "TROOT.h"
#include "Index.hpp"

GenPartIndex IsDiHiggsSL(Int_t const* pdgIds, Int_t const* motherIdxs, Int_t const* statuses, Int_t nGenPart);
GenJetIndex Match(GenPartIndex const& genPartIdx, PtEtaPhiMArray const& genPart, PtEtaPhiMArray const& genJet, Int_t const* genPartId, Int_t const* jetPartFlav);

#endif