#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <vector>
#include <map>

#ifndef N_LEP
#define N_LEP 1
#endif

inline constexpr size_t MAX_GEN_JET = 20;
inline constexpr size_t MAX_RECO_JET = 12;
inline constexpr size_t MAX_RECO_LEP = 2;
inline constexpr size_t MAX_GEN_LEP = 2;

enum class LepIdx { lep1, lep2, count };

// objects in SL channel resolved topology
enum class ObjSLRes { b1, b2, q1, q2, lep, met, count };

// objects in DL channel resolved topology
enum class ObjDLRes { b1, b2, lep1, lep2, met, count };

// used to index 4-vectors in the input of HME and 4-vectors in gen truth
using Obj = std::conditional_t<N_LEP == 1, ObjSLRes, ObjDLRes>;

// mapping Obj to corresponding gen truth branch names
// static const std::map<ObjSLRes, char const*> GenTruthMapSL = { { ObjSLRes::b1, "genb1" },
//                                                                { ObjSLRes::b2, "genb2" },
//                                                                { ObjSLRes::q1, "genV2prod1" },
//                                                                { ObjSLRes::q2, "genV2prod2" },
//                                                                { ObjSLRes::lep, "genV1prod1" } };

// static const std::map<ObjDLRes, char const*> GenTruthMapDL = { { ObjDLRes::b1, "genb1" },
//                                                                { ObjDLRes::b2, "genb2" },
//                                                                { ObjDLRes::lep1, "genV2prod1" },
//                                                                { ObjDLRes::lep2, "genV2prod1" } };

// static const std::map<Obj, char const*> GenTruthMap = N_LEP == 1 ? GenTruthMapSL : GenTruthMapDL;

// PDFs in SL channel resolved topology
enum class PDF1dSLRes { pdf_b1, count };
enum class PDF2dSLRes { pdf_b1b2, count };
static const std::vector<char const*> pdf1d_names { "pdf_b1" };
static const std::vector<char const*> pdf2d_names { "pdf_b1b2" };

#endif