#ifndef READOUT_CONST_HPP
#define READOUT_CONST_HPP

namespace ReadConst
{
    enum class Lep { lep1, lep2, count };
    enum class Nu { nu1, nu2, count };
    enum class Quark { b1, b2, q1, q2, count };
    inline constexpr size_t MAX_GEN_LEP = static_cast<size_t>(Lep::count);
    inline constexpr size_t MAX_GEN_QUARK = static_cast<size_t>(Quark::count);
    inline constexpr size_t MAX_GEN_NU = static_cast<size_t>(Nu::count);
    
    // inline constexpr size_t MAX_GEN_JET = 20;
    inline constexpr size_t MAX_RECO_JET = 20;
    inline constexpr size_t MAX_RECO_FAT_JET = 5;
    inline constexpr size_t MAX_RECO_LEP = 2;
    inline constexpr size_t NUM_BQ = 2;
    inline constexpr size_t NUM_LQ = 2;
}; 

#endif