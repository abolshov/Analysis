#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <vector>
#include <map>

#ifndef N_LEP
#define N_LEP 1
#endif

enum class Channel { SL, DL };

inline constexpr size_t MAX_GEN_JET = 20;
inline constexpr size_t MAX_RECO_JET = 12;
inline constexpr size_t MAX_RECO_LEP = 2;
inline constexpr size_t MAX_GEN_LEP = 2;

// obj name to its index in HME input (or index in gen truth)
static const std::map<std::string, size_t> GenTruthIdxMapSL = { { "b1", 0 },
                                                                { "b2", 1 },
                                                                { "q1", 2 },
                                                                { "q2", 3 },
                                                                { "lep", 4 },
                                                                { "met", 5 } };

static const std::map<std::string, size_t> GenTruthIdxMapDL = { { "b1", 0 },
                                                                { "b2", 1 },
                                                                { "lep1", 2 },
                                                                { "lep2", 3 },
                                                                { "met", 4 } };

// obj name to name of the branch with pt, eta, phi, mass of that object
static const std::map<std::string, std::string> GenTruthBranchMapSL = { { "b1", "genb1" },
                                                                        { "b2", "genb2" },
                                                                        { "q1", "genV2prod1" },
                                                                        { "q2", "genV2prod2" },
                                                                        { "lep", "genV1prod1" },
                                                                        { "met", "GenMET" } };

static const std::map<std::string, std::string> GenTruthBranchMapDL = { { "b1", "genb1" },
                                                                        { "b2", "genb2" },
                                                                        { "lep1", "genV1prod1" },
                                                                        { "lep2", "genV2prod1" },
                                                                        { "met", "GenMET" } };

static const std::vector<std::string> KinVarNames = { "_pt", "_eta", "_phi", "_mass" };

// PDFs in SL channel resolved topology
enum class PDF1dSLRes { pdf_b1, count };
enum class PDF2dSLRes { pdf_b1b2, count };
static const std::vector<char const*> pdf1d_names { "pdf_b1" };
static const std::vector<char const*> pdf2d_names { "pdf_b1b2" };

#endif