#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>

enum class Channel { SL, DL };
enum class Topology { Resolved, Boosted };
enum class Mode { Validation, Estimation };

enum class Lep { lep1, lep2 };
enum class Nu { nu1, nu2 };
enum class Quark { b1, b2, q1, q2 };

inline constexpr size_t MAX_GEN_JET = 20;
inline constexpr size_t MAX_RECO_JET = 12;
inline constexpr size_t MAX_RECO_LEP = 2;
inline constexpr size_t MAX_GEN_LEP = 2;
inline constexpr size_t MAX_GEN_QUARK = 4;
inline constexpr size_t MAX_GEN_NU = 2;

inline constexpr double HIGGS_MASS = 125.03;
inline constexpr double HIGGS_WIDTH = 0.004;
inline constexpr double TOL = 10e-7;
inline constexpr int N_ATTEMPTS = 1;
inline constexpr int N_ITER = 1000;

inline constexpr double MAX_MASS = 2000;
inline constexpr int N_BINS = 5000;

inline constexpr double MET_SIGMA = 25.2;

// PDFs in SL channel resolved topology
enum class PDF1dSLRes { pdf_b1, count };
enum class PDF2dSLRes { pdf_b1b2, count };

#endif