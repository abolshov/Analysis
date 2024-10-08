#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include <optional>
#include <vector>
#include <memory>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

// general parameters
inline constexpr double HIGGS_MASS = 125.03;
inline constexpr double HIGGS_WIDTH = 0.004;
inline constexpr double TOL = 10e-7;
inline constexpr int N_ATTEMPTS = 1;
inline constexpr int N_ITER = 1000;

// event mass distribution parameters
inline constexpr double MAX_MASS = 2000;
inline constexpr int N_BINS = 5000;

inline constexpr double MET_SIGMA = 25.2;

using OptionalPair = std::optional<std::pair<double, double>>;

enum PhysObj { bj1, bj2, wj1, wj2, lep, met };

// computes rescaling factors for jets using rescale pdf and mass of particle->jj
// has side effect: may swap j1 and j2
OptionalPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass);

// compute 4-momentum of neutrino from W->lv using Higgs H->WW mass constraint
std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, double mh, double eta);

// compute 4-momentum of neutrino from W->lv using Higgs W->lv mass constraint
std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& met, double mw, double numet_dphi);

// this is HME: samples pdfs, computes corrections and calculates mass distribution for an event
// in case of success returns mass and fraction of iterations which computed mass
OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt, std::pair<double, double> lj_pt_res);

#endif
