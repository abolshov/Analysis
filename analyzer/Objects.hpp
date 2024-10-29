#ifndef OBJECTS_HPP
#define OBJECTS_HPP

#include <memory>
#include <iostream>

#include "TROOT.h"

#include "Constants.hpp"

using UArrF_t = std::unique_ptr<Float_t[]>;
using UArrUC_t = std::unique_ptr<UChar_t[]>;
using UArrS_t = std::unique_ptr<Short_t[]>;
using UArrI_t = std::unique_ptr<Int_t[]>;

struct Kinematics_t
{
	Kinematics_t(size_t sz);

	UArrF_t pt;
	UArrF_t eta;
	UArrF_t phi;
	UArrF_t mass;
};

struct GenJet_t : public Kinematics_t 
{
    GenJet_t();

    Int_t nGenJet;
    UArrS_t part_flav;
    UArrUC_t hadr_flav;
};

struct RecoJet_t : public Kinematics_t 
{
    RecoJet_t();

    Int_t nRecoJet;
    UArrUC_t part_flav;
    UArrS_t hadr_flav;

    UArrF_t btagPNetB;
    UArrF_t PNetRegPtRawCorr;
    UArrF_t PNetRegPtRawRes;
    UArrF_t PNetRegPtRawCorrNu;
};

struct RecoLep_t : public Kinematics_t
{
    RecoLep_t();

    UArrF_t lep_iso;
    UArrI_t lep_type;
};

struct GenLep_t : public Kinematics_t
{
    GenLep_t();

    UArrI_t pdgId;
};

struct Particle 
{
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t mass;
};

#endif