#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// Event selection constants
inline constexpr int MAX_GENPART = 250;
inline constexpr int MAX_AK4_GENJET = 20;
inline constexpr int N_SIG_PART = 16;
inline constexpr int N_B_JETS = 2;
inline constexpr int N_LIGHT_JETS = 2;

inline constexpr int RADION_ID = 35;

enum Signal { Xfirst, Xlast, HbbFirst, HWWfirst, HbbLast, HWWlast, b, bbar, LepWfirst, HadWfirst, HadWlast, q1, q2, LepWlast, l, nu };

#endif