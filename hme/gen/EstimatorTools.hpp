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
inline constexpr double MAX_MASS = 1500;
inline constexpr int N_BINS = 1500;

using OptionalPair = std::optional<std::pair<double, double>>;

enum PhysObj { bj1, bj2, wj1, wj2, lep, met };
enum PDF { lead_bjet_pt_pdf, subl_bjet_pt_pdf, lead_bjet_E_pdf, subl_bjet_E_pdf, lead_bjet_pz_pdf, subl_bjet_pz_pdf };

// computes rescaling factors for jets using rescale pdf and mass of particle->jj
// has side effect: may swap j1 and j2
OptionalPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass);

// compute 4-momentum of neutrino from W->lv using Higgs H->WW mass constraint
std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, double mh, double eta);
std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, double mh, double eta, double phi);

// this is HME: samples pdfs, computes corrections and calculates mass distribution for an event
// in case of success returns mass and fraction of iterations which computed mass
OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt);

namespace Experimental
{
    OptionalPair EstimateMassIdealEtaPhi(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt, std::pair<double, double> const& dir);
    OptionalPair EstimateMassIdealNu(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& lead_b_pdf, TRandom3& rg, int evt, TLorentzVector const& nu);
    OptionalPair EstimateMassIdealNu2dPDF(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH2F>& bjet_2d_pdf, TRandom3& rg, int evt, TLorentzVector const& nu);
    OptionalPair EstimateMassNew(std::vector<TLorentzVector> const& particles, std::vector<std::unique_ptr<TH1F>> const& pdfs, TRandom3& rg, int evt, TLorentzVector const& nu);
    OptionalPair EstimateMassIdealHbb(std::vector<TLorentzVector> const& particles, TRandom3& rg, int evt); 
    OptionalPair EstimateMassIdealHWW(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& pdf, TRandom3& rg, int evt, TLorentzVector const& HtoWW);    
}

#endif