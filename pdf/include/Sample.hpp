#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "TString.h"
#include "Constants.hpp"

struct Sample
{
    TString file_name{};
    TString tree_name{"Events"};
    TString type;
    Channel channel;
    Int_t masspoint = -1;
};

#endif