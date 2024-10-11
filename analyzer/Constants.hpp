#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <vector>

// namespace for gen constants 
namespace Gen 
{
    inline constexpr size_t MAX_GEN_JET = 20;
};

// namespace for reco constants 
namespace Reco 
{
    inline constexpr size_t MAX_RECO_JET = 12;
};

// objects in SL channel resolved topology
enum ObjSLRes { b1, b2, q1, q2, lep, met, count };

// PDFs in SL channel resolved topology
enum PDF1dSLRes { pdf_b1, count };
enum PDF2dSLRes { pdf_b1b2, count };
static const std::vector<char const*> pdf1d_names { "pdf_b1" };
static const std::vector<char const*> pdf2d_names { "pdf_b1b2" };

#endif