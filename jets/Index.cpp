#include "Index.hpp"

std::ostream& operator<<(std::ostream& os, GenPartIndex const& idx)
{
    auto const& [h1, h2, w1, w2, b1, b2, q1, q2, l, nu] = idx;
    os << "(" << h1 << ", "
       << h2 << ", "
       << w1 << ", "
       << w2 << ", "
       << b1 << ", "
       << b2 << ", "
       << q1 << ", "
       << q2 << ", "
       << l << ", "
       << nu << ")"; 
    return os;
}

std::ostream& operator<<(std::ostream& os, GenJetIndex const& idx)
{
    auto const& [bj1, bj2, lj1, lj2] = idx;
    os << "(" << bj1 << ", "
       << bj2 << ", "
       << lj1 << ", "
       << lj2 << ")"; 
    return os;
}