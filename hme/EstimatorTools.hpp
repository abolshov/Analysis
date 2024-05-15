#ifndef ESTIMATOR_TOOLS_HPP
#define ESTIMATOR_TOOLS_HPP

#include <optional>
#include <vector>
#include <memory>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TH1.h"

// general parameters
static constexpr double HIGGS_MASS = 125.03;
static constexpr double HIGGS_WIDTH = 0.004;
static constexpr double TOL = 10e-7;
static constexpr int N_ATTEMPTS = 100;
static constexpr int N_ITER = 1000;

// event mass distribution parameters
static constexpr double MAX_MASS = 1500;
static constexpr int N_BINS = 1500;

using OptionalPair = std::optional<std::pair<double, double>>;

enum PhysObj { bj1, bj2, wj1, wj2, lep, met };

// equation for pz of neutrino is quadratic, so 2 solutions are possible
// if discriminant < 0 returns nullopt
// otherwise returns pair of solutions
OptionalPair AnalyticalMass(std::vector<TLorentzVector> const& particles);

// computes rescaling factors for jets using rescale pdf and mass of particle->jj
// has side effect: may swap j1 and j2
OptionalPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass);

// compute 4-momentum of neutrino from W->lv using Higgs H->WW mass constraint
std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, double mh, double eta);

// this is HME: samples pdfs, computes corrections and calculates mass distribution for an event
// in case of success returns mass and fraction of iterations which computed mass
OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt);

#endif