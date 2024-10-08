#ifndef OBJECTS_HPP
#define OBJECTS_HPP

#include <memory>

#include "TROOT.h"

#include "Constants.hpp"

using UArrF_t = std::unique_ptr<Float_t[]>;
using UArrUC_t = std::unique_ptr<UChar_t[]>;
using UArrS_t = std::unique_ptr<Short_t[]>;

struct Kinematics
{
	Kinematics(size_t sz);

	UArrF_t pt;
	UArrF_t eta;
	UArrF_t phi;
	UArrF_t mass;
};

struct GenJet : public Kinematics 
{
    GenJet();

    Int_t nGenJet;
    // UArrUC_t part_flav;
    // UArrS_t hadr_flav;
};

struct RecoJet : public Kinematics 
{
    RecoJet();

    Int_t nRecoJet;
    UArrUC_t part_flav;
    UArrS_t hadr_flav;

    UArrF_t btagPNetB;
    UArrF_t PNetRegPtRawCorr;
    UArrF_t PNetRegPtRawRes;
    UArrF_t PNetRegPtRawCorrNu;
};

#endif