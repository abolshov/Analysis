#include "Objects.hpp"

Kinematics_t::Kinematics_t(size_t sz)
:	pt(std::make_unique<Float_t[]>(sz))
,	eta(std::make_unique<Float_t[]>(sz))
,	phi(std::make_unique<Float_t[]>(sz))
,	mass(std::make_unique<Float_t[]>(sz))
{
    for (size_t i = 0; i < sz; ++i)
	{
		pt[i] = 0.0;
		eta[i] = 0.0;
		phi[i] = 0.0;
		mass[i] = 0.0;
	}
}

GenJet_t::GenJet_t() 
:	Kinematics_t(MAX_GEN_JET)
,	nGenJet(0)
,	part_flav(std::make_unique<Short_t[]>(MAX_GEN_JET))
,	hadr_flav(std::make_unique<UChar_t[]>(MAX_GEN_JET))
{
	for (size_t i = 0; i < MAX_GEN_JET; ++i)
	{
		part_flav[i] = 0;
		hadr_flav[i] = 0;
	}
}

RecoJet_t::RecoJet_t() 
: 	Kinematics_t(MAX_RECO_JET)
,	nRecoJet(0)
,	part_flav(std::make_unique<UChar_t[]>(MAX_RECO_JET))
,	hadr_flav(std::make_unique<Short_t[]>(MAX_RECO_JET))
,	btagPNetB(std::make_unique<Float_t[]>(MAX_RECO_JET))
,	PNetRegPtRawCorr(std::make_unique<Float_t[]>(MAX_RECO_JET))
,	PNetRegPtRawRes(std::make_unique<Float_t[]>(MAX_RECO_JET))
,	PNetRegPtRawCorrNu(std::make_unique<Float_t[]>(MAX_RECO_JET))
{
	for (size_t i = 0; i < MAX_RECO_JET; ++i)
	{
		part_flav[i] = static_cast<UChar_t>(0);
		hadr_flav[i] = 0;
		btagPNetB[i] = 0.0;
		PNetRegPtRawCorr[i] = 0.0;
		PNetRegPtRawRes[i] = 0.0;
		PNetRegPtRawCorrNu[i] = 0.0;
	}
}

RecoLep_t::RecoLep_t() 
:	Kinematics_t(MAX_RECO_LEP)
,	lep_iso(std::make_unique<Float_t[]>(MAX_RECO_LEP))
,	lep_type(std::make_unique<Int_t[]>(MAX_RECO_LEP))
{
	for (size_t i = 0; i < MAX_RECO_LEP; ++i)
	{
		lep_type[i] = 0;
		lep_iso[i] = 0.0;
	}
}

GenLep_t::GenLep_t() 
:	Kinematics_t(MAX_GEN_LEP)
,	pdgId(std::make_unique<Int_t[]>(MAX_GEN_LEP))
{
	for (size_t i = 0; i < MAX_GEN_LEP; ++i)
	{
		pdgId[i] = 0;
	}
}