#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>

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

#endif