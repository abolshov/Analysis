#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <vector>

inline constexpr size_t MAX_GEN_JET = 20;
inline constexpr size_t MAX_RECO_JET = 12;
inline constexpr size_t MAX_RECO_LEP = 2;
inline constexpr size_t MAX_GEN_LEP = 2;

enum class LepIdx { lep1, lep2, count };

// objects in SL channel resolved topology
enum class ObjSLRes { b1, b2, q1, q2, lep, met, count };

// PDFs in SL channel resolved topology
enum class PDF1dSLRes { pdf_b1, count };
enum class PDF2dSLRes { pdf_b1b2, count };
static const std::vector<char const*> pdf1d_names { "pdf_b1" };
static const std::vector<char const*> pdf2d_names { "pdf_b1b2" };

#endif