#ifndef JET_COMB_HPP
#define JET_COMB_HPP

#include "TString.h"

struct JetComb
{
    int b1 = -1;
    int b2 = -1;
    int q1 = -1;
    int q2 = -1;

    inline TString ToTString()
    {
        return q1 < 0 ? Form("b%db%d", b1, b2) : Form("b%db%dq%dq%d", b1, b2, q1, q2);
    }
};

#endif