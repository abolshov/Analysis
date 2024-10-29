#ifndef DEF_HPP
#define DEF_HPP

#include <vector>
#include <memory>

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"

#include "TH1.h"
#include "TH2.h"

using UHist1D = std::unique_ptr<TH1F>;
using UHist2D = std::unique_ptr<TH2F>;

using HistVec1d_t = std::vector<UHist1D>;
using HistVec2d_t = std::vector<UHist2D>;

using LorentzVectorF_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Float_t>>;
using VecLVF_t = std::vector<LorentzVectorF_t>;

#endif