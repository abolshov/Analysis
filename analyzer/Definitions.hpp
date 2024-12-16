#ifndef DEF_HPP
#define DEF_HPP

#include <vector>
#include <memory>
#include <optional>

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"

#include "TH1.h"
#include "TH2.h"

template <typename T> 
using UHist_t = std::unique_ptr<T>;

template <typename T> 
using HistVec_t = std::vector<UHist_t<T>>;

using LorentzVectorF_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Float_t>>;
using VecLVF_t = std::vector<LorentzVectorF_t>;

using OptionalPair = std::optional<std::pair<double, double>>;

#endif