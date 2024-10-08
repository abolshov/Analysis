#ifndef OBJECTS_HPP
#define OBJECTS_HPP

#include <memory>

#include "TROOT.h"
#include "TLorentzVector.h"

#include "Constants.hpp"

using UArrF_t = std::unique_ptr<Float_t[]>;
using UArrUC_t = std::unique_ptr<UChar_t[]>;
using UArrS_t = std::unique_ptr<Short_t[]>;
using UArrLV_t = std::unique_ptr<TLorentzVector[]>;

struct Kinematics
{
	Kinematics(size_t sz);

	UArrF_t pt;
	UArrF_t eta;
	UArrF_t phi;
	UArrF_t mass;
    UArrLV_t p4;

    protected:
    void SetP4(size_t sz);
};

struct GenJet : public Kinematics 
{
    GenJet();
    inline void SetP4() { Kinematics::SetP4(nGenJet); }

    Int_t nGenJet;
    // UArrUC_t part_flav;
    // UArrS_t hadr_flav;
};

struct RecoJet : public Kinematics 
{
    RecoJet();
    inline void SetP4() { Kinematics::SetP4(nRecoJet); }

    Int_t nRecoJet;
    UArrUC_t part_flav;
    UArrS_t hadr_flav;

    UArrF_t btagPNetB;
    UArrF_t PNetRegPtRawCorr;
    UArrF_t PNetRegPtRawRes;
    UArrF_t PNetRegPtRawCorrNu;
};

#endif