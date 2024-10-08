#include "Objects.hpp"

Kinematics::Kinematics(size_t sz) 
: pt(std::make_unique<Float_t[]>(sz)),
  eta(std::make_unique<Float_t[]>(sz)),
  phi(std::make_unique<Float_t[]>(sz)),
  mass(std::make_unique<Float_t[]>(sz))
{}

GenJet::GenJet() 
: Kinematics(Gen::MAX_GEN_JET), 
  nGenJet(0)
//   part_flav(std::make_unique<UChar_t>(Gen::MAX_GEN_JET)),
//   hadr_flav(std::make_unique<Short_t>(Gen::MAX_GEN_JET))
{}

RecoJet::RecoJet() 
: Kinematics(Reco::MAX_RECO_JET),
  nRecoJet(0),
  part_flav(std::make_unique<UChar_t[]>(Reco::MAX_RECO_JET)),
  hadr_flav(std::make_unique<Short_t[]>(Reco::MAX_RECO_JET)),
  btagPNetB(std::make_unique<Float_t[]>(Reco::MAX_RECO_JET)),
  PNetRegPtRawCorr(std::make_unique<Float_t[]>(Reco::MAX_RECO_JET)),
  PNetRegPtRawRes(std::make_unique<Float_t[]>(Reco::MAX_RECO_JET)),
  PNetRegPtRawCorrNu(std::make_unique<Float_t[]>(Reco::MAX_RECO_JET))
{}
