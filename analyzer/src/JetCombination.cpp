#include "JetCombination.hpp"

#ifdef DEV
    #include <unordered_set>

    bool JetComb::HasUniqueJets(Channel ch) const
    {
        std::unordered_set<int> indices;
        indices.insert(b1);
        indices.insert(b2);
        if (ch == Channel::SL)
        {
            indices.insert(q1);
            indices.insert(q2);
        }
        // return indices.size() == QuarkCount(ch) && !indices.contains(-1); // c++20
        return indices.size() == QuarkCount(ch) && !indices.count(-1);
    }
#endif