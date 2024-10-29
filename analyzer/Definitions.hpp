#ifndef DEF_HPP
#define DEF_HPP

#include <vector>
#include <memory>
#include <optional>

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"

#include "TH1.h"
#include "TH2.h"

using UHist1d_t = std::unique_ptr<TH1F>;
using UHist2d_t = std::unique_ptr<TH2F>;

using HistVec1d_t = std::vector<UHist1d_t>;
using HistVec2d_t = std::vector<UHist2d_t>;

using LorentzVectorF_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Float_t>>;
using VecLVF_t = std::vector<LorentzVectorF_t>;

using OptionalPair = std::optional<std::pair<double, double>>;

#endif