#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>

inline constexpr int SEED = 42;

enum class Channel { SL, DL };
enum class Topology { Resolved, Boosted };
enum class Mode { Validation, Estimation };

enum class Lep { lep1, lep2, count };
enum class Nu { nu1, nu2, count };
enum class Quark { b1, b2, q1, q2, count };
inline constexpr size_t MAX_GEN_LEP = static_cast<size_t>(Lep::count);
inline constexpr size_t MAX_GEN_QUARK = static_cast<size_t>(Quark::count);
inline constexpr size_t MAX_GEN_NU = static_cast<size_t>(Nu::count);

// PDFs in SL channel resolved topology
enum class PDF1 { numet_pt, numet_dphi, nulep_deta, hh_dphi, mbb, mww, hh_deta, count };
enum class PDF2 { b1b2, hh_dEtadPhi, hh_pt_e, count };
inline constexpr size_t NUM_PDF_1D = static_cast<size_t>(PDF1::count);
inline constexpr size_t NUM_PDF_2D = static_cast<size_t>(PDF2::count);

// return values of Estimator
enum class Output { mass, integral, width, peak_val, count };
inline constexpr size_t OUTPUT_SIZE = static_cast<size_t>(Output::count);

// objects
enum class ObjSL { bj1, bj2, lj1, lj2, lep, met, count };

inline constexpr size_t MAX_GEN_JET = 20;
inline constexpr size_t MAX_RECO_JET = 12;
inline constexpr size_t MAX_RECO_LEP = 2;
inline constexpr size_t NUM_BQ = 2;
inline constexpr size_t NUM_LQ = 2;

inline constexpr Float_t HIGGS_MASS = 125.03;
inline constexpr Float_t HIGGS_WIDTH = 0.004;
inline constexpr Float_t TOL = 10e-7;
inline constexpr int N_ATTEMPTS = 1;
inline constexpr int N_ITER = 1000;

inline constexpr Float_t MAX_MASS = 5000.0;
inline constexpr Float_t MIN_MASS = 200.0;
inline constexpr int N_BINS = 10000;

inline constexpr Float_t MET_SIGMA = 25.2;

inline constexpr unsigned Q16 = 16;
inline constexpr unsigned Q84 = 84;

inline constexpr size_t NUM_BEST_BTAG = 3;

#endif