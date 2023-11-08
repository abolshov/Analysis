//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 24 09:41:28 2023 by ROOT version 6.28/06
// from TTree syncTree/Friend tree for Events
// found on file: NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root
//////////////////////////////////////////////////////////

#ifndef syncTree_h
#define syncTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class syncTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         puWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   UInt_t          nJet;
   Float_t         Jet_btagSF_deepjet_M_down[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_M[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_M_up[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_hf[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_cferr1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_jes[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_cferr2[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_lf[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_lf[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_cferr1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_lfstats1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_lfstats2[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_hfstats1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_hfstats2[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_lfstats2[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_hf[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_lfstats1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_jes[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_hfstats2[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_down_hfstats1[25];   //[nJet]
   Float_t         Jet_btagSF_deepjet_shape_up_cferr2[25];   //[nJet]
   UInt_t          nGenHiggs;
   Int_t           GenHiggs_status[2];   //[nGenHiggs]
   Float_t         GenHiggs_phi[2];   //[nGenHiggs]
   Float_t         GenHiggs_eta[2];   //[nGenHiggs]
   Float_t         GenHiggs_mass[2];   //[nGenHiggs]
   Int_t           GenHiggs_idx[2];   //[nGenHiggs]
   Float_t         GenHiggs_pt[2];   //[nGenHiggs]
   Int_t           GenHiggs_statusFlags[2];   //[nGenHiggs]
   Int_t           GenHiggs_pdgId[2];   //[nGenHiggs]
   Int_t           GenHiggs_charge[2];   //[nGenHiggs]
   UInt_t          nGenLepFromW1FromHiggs;
   Int_t           GenLepFromW1FromHiggs_status[1];   //[nGenLepFromW1FromHiggs]
   Float_t         GenLepFromW1FromHiggs_phi[1];   //[nGenLepFromW1FromHiggs]
   Float_t         GenLepFromW1FromHiggs_eta[1];   //[nGenLepFromW1FromHiggs]
   Float_t         GenLepFromW1FromHiggs_mass[1];   //[nGenLepFromW1FromHiggs]
   Int_t           GenLepFromW1FromHiggs_idx[1];   //[nGenLepFromW1FromHiggs]
   Float_t         GenLepFromW1FromHiggs_pt[1];   //[nGenLepFromW1FromHiggs]
   Int_t           GenLepFromW1FromHiggs_statusFlags[1];   //[nGenLepFromW1FromHiggs]
   Int_t           GenLepFromW1FromHiggs_pdgId[1];   //[nGenLepFromW1FromHiggs]
   Int_t           GenLepFromW1FromHiggs_charge[1];   //[nGenLepFromW1FromHiggs]
   UInt_t          nGenX;
   Int_t           GenX_status[1];   //[nGenX]
   Float_t         GenX_phi[1];   //[nGenX]
   Float_t         GenX_eta[1];   //[nGenX]
   Float_t         GenX_mass[1];   //[nGenX]
   Int_t           GenX_idx[1];   //[nGenX]
   Float_t         GenX_pt[1];   //[nGenX]
   Int_t           GenX_statusFlags[1];   //[nGenX]
   Int_t           GenX_pdgId[1];   //[nGenX]
   Int_t           GenX_charge[1];   //[nGenX]
   UInt_t          nGenLepFromW2FromHiggs;
   Int_t           GenLepFromW2FromHiggs_status[1];   //[nGenLepFromW2FromHiggs]
   Float_t         GenLepFromW2FromHiggs_phi[1];   //[nGenLepFromW2FromHiggs]
   Float_t         GenLepFromW2FromHiggs_eta[1];   //[nGenLepFromW2FromHiggs]
   Float_t         GenLepFromW2FromHiggs_mass[1];   //[nGenLepFromW2FromHiggs]
   Int_t           GenLepFromW2FromHiggs_idx[1];   //[nGenLepFromW2FromHiggs]
   Float_t         GenLepFromW2FromHiggs_pt[1];   //[nGenLepFromW2FromHiggs]
   Int_t           GenLepFromW2FromHiggs_statusFlags[1];   //[nGenLepFromW2FromHiggs]
   Int_t           GenLepFromW2FromHiggs_pdgId[1];   //[nGenLepFromW2FromHiggs]
   Int_t           GenLepFromW2FromHiggs_charge[1];   //[nGenLepFromW2FromHiggs]
   UInt_t          nGenW1FromHiggs;
   Int_t           GenW1FromHiggs_status[1];   //[nGenW1FromHiggs]
   Float_t         GenW1FromHiggs_phi[1];   //[nGenW1FromHiggs]
   Float_t         GenW1FromHiggs_eta[1];   //[nGenW1FromHiggs]
   Float_t         GenW1FromHiggs_mass[1];   //[nGenW1FromHiggs]
   Int_t           GenW1FromHiggs_idx[1];   //[nGenW1FromHiggs]
   Float_t         GenW1FromHiggs_pt[1];   //[nGenW1FromHiggs]
   Int_t           GenW1FromHiggs_statusFlags[1];   //[nGenW1FromHiggs]
   Int_t           GenW1FromHiggs_pdgId[1];   //[nGenW1FromHiggs]
   Int_t           GenW1FromHiggs_charge[1];   //[nGenW1FromHiggs]
   UInt_t          nGenBQuarkFromHiggs;
   Int_t           GenBQuarkFromHiggs_status[2];   //[nGenBQuarkFromHiggs]
   Float_t         GenBQuarkFromHiggs_phi[2];   //[nGenBQuarkFromHiggs]
   Float_t         GenBQuarkFromHiggs_eta[2];   //[nGenBQuarkFromHiggs]
   Float_t         GenBQuarkFromHiggs_mass[2];   //[nGenBQuarkFromHiggs]
   Int_t           GenBQuarkFromHiggs_idx[2];   //[nGenBQuarkFromHiggs]
   Float_t         GenBQuarkFromHiggs_pt[2];   //[nGenBQuarkFromHiggs]
   Int_t           GenBQuarkFromHiggs_statusFlags[2];   //[nGenBQuarkFromHiggs]
   Int_t           GenBQuarkFromHiggs_pdgId[2];   //[nGenBQuarkFromHiggs]
   Int_t           GenBQuarkFromHiggs_charge[2];   //[nGenBQuarkFromHiggs]
   UInt_t          nGenW2FromHiggs;
   Int_t           GenW2FromHiggs_status[1];   //[nGenW2FromHiggs]
   Float_t         GenW2FromHiggs_phi[1];   //[nGenW2FromHiggs]
   Float_t         GenW2FromHiggs_eta[1];   //[nGenW2FromHiggs]
   Float_t         GenW2FromHiggs_mass[1];   //[nGenW2FromHiggs]
   Int_t           GenW2FromHiggs_idx[1];   //[nGenW2FromHiggs]
   Float_t         GenW2FromHiggs_pt[1];   //[nGenW2FromHiggs]
   Int_t           GenW2FromHiggs_statusFlags[1];   //[nGenW2FromHiggs]
   Int_t           GenW2FromHiggs_pdgId[1];   //[nGenW2FromHiggs]
   Int_t           GenW2FromHiggs_charge[1];   //[nGenW2FromHiggs]
   UInt_t          nGenNuFromW2FromHiggs;
   Int_t           GenNuFromW2FromHiggs_status[1];   //[nGenNuFromW2FromHiggs]
   Float_t         GenNuFromW2FromHiggs_phi[1];   //[nGenNuFromW2FromHiggs]
   Float_t         GenNuFromW2FromHiggs_eta[1];   //[nGenNuFromW2FromHiggs]
   Float_t         GenNuFromW2FromHiggs_mass[1];   //[nGenNuFromW2FromHiggs]
   Int_t           GenNuFromW2FromHiggs_idx[1];   //[nGenNuFromW2FromHiggs]
   Float_t         GenNuFromW2FromHiggs_pt[1];   //[nGenNuFromW2FromHiggs]
   Int_t           GenNuFromW2FromHiggs_statusFlags[1];   //[nGenNuFromW2FromHiggs]
   Int_t           GenNuFromW2FromHiggs_pdgId[1];   //[nGenNuFromW2FromHiggs]
   Int_t           GenNuFromW2FromHiggs_charge[1];   //[nGenNuFromW2FromHiggs]
   UInt_t          nGenQuarkFromW2FromHiggs;
   Int_t           GenQuarkFromW2FromHiggs_status[2];   //[nGenQuarkFromW2FromHiggs]
   Float_t         GenQuarkFromW2FromHiggs_phi[2];   //[nGenQuarkFromW2FromHiggs]
   Float_t         GenQuarkFromW2FromHiggs_eta[2];   //[nGenQuarkFromW2FromHiggs]
   Float_t         GenQuarkFromW2FromHiggs_mass[2];   //[nGenQuarkFromW2FromHiggs]
   Int_t           GenQuarkFromW2FromHiggs_idx[2];   //[nGenQuarkFromW2FromHiggs]
   Float_t         GenQuarkFromW2FromHiggs_pt[2];   //[nGenQuarkFromW2FromHiggs]
   Int_t           GenQuarkFromW2FromHiggs_statusFlags[2];   //[nGenQuarkFromW2FromHiggs]
   Int_t           GenQuarkFromW2FromHiggs_pdgId[2];   //[nGenQuarkFromW2FromHiggs]
   Int_t           GenQuarkFromW2FromHiggs_charge[2];   //[nGenQuarkFromW2FromHiggs]
   UInt_t          nGenNuFromW1FromHiggs;
   Int_t           GenNuFromW1FromHiggs_status[1];   //[nGenNuFromW1FromHiggs]
   Float_t         GenNuFromW1FromHiggs_phi[1];   //[nGenNuFromW1FromHiggs]
   Float_t         GenNuFromW1FromHiggs_eta[1];   //[nGenNuFromW1FromHiggs]
   Float_t         GenNuFromW1FromHiggs_mass[1];   //[nGenNuFromW1FromHiggs]
   Int_t           GenNuFromW1FromHiggs_idx[1];   //[nGenNuFromW1FromHiggs]
   Float_t         GenNuFromW1FromHiggs_pt[1];   //[nGenNuFromW1FromHiggs]
   Int_t           GenNuFromW1FromHiggs_statusFlags[1];   //[nGenNuFromW1FromHiggs]
   Int_t           GenNuFromW1FromHiggs_pdgId[1];   //[nGenNuFromW1FromHiggs]
   Int_t           GenNuFromW1FromHiggs_charge[1];   //[nGenNuFromW1FromHiggs]
   Int_t           event;
   Int_t           ls;
   Int_t           run;
   Int_t           n_presel_mu;
   Int_t           n_fakeablesel_mu;
   Int_t           n_mvasel_mu;
   Int_t           n_presel_ele;
   Int_t           n_fakeablesel_ele;
   Int_t           n_mvasel_ele;
   Int_t           n_presel_ak4Jet;
   Int_t           n_presel_ak4JetVBF;
   Int_t           n_presel_ak8Jet;
   Int_t           n_presel_ak8lsJet;
   Int_t           n_loose_ak4BJet;
   Int_t           n_medium_ak4BJet;
   Int_t           is_ee;
   Int_t           is_mm;
   Int_t           is_em;
   Int_t           is_boosted;
   Int_t           is_semiboosted;
   Int_t           is_resolved;
   Float_t         mu1_pt;
   Float_t         mu1_conept;
   Float_t         mu1_eta;
   Float_t         mu1_phi;
   Float_t         mu1_E;
   Int_t           mu1_charge;
   Float_t         mu1_miniRelIso;
   Float_t         mu1_PFRelIso04;
   Int_t           mu1_jetNDauChargedMVASel;
   Float_t         mu1_jetPtRel;
   Float_t         mu1_jetRelIso;
   Float_t         mu1_jetDeepJet;
   Float_t         mu1_sip3D;
   Float_t         mu1_dxy;
   Float_t         mu1_dxyAbs;
   Float_t         mu1_dz;
   Float_t         mu1_segmentCompatibility;
   Float_t         mu1_leptonMVA;
   Int_t           mu1_mediumID;
   Float_t         mu1_dpt_div_pt;
   Int_t           mu1_isfakeablesel;
   Int_t           mu1_ismvasel;
   Float_t         mu1_isGenMatched;
   Float_t         mu2_pt;
   Float_t         mu2_conept;
   Float_t         mu2_eta;
   Float_t         mu2_phi;
   Float_t         mu2_E;
   Int_t           mu2_charge;
   Float_t         mu2_miniRelIso;
   Float_t         mu2_PFRelIso04;
   Int_t           mu2_jetNDauChargedMVASel;
   Float_t         mu2_jetPtRel;
   Float_t         mu2_jetRelIso;
   Float_t         mu2_jetDeepJet;
   Float_t         mu2_sip3D;
   Float_t         mu2_dxy;
   Float_t         mu2_dxyAbs;
   Float_t         mu2_dz;
   Float_t         mu2_segmentCompatibility;
   Float_t         mu2_leptonMVA;
   Int_t           mu2_mediumID;
   Float_t         mu2_dpt_div_pt;
   Int_t           mu2_isfakeablesel;
   Int_t           mu2_ismvasel;
   Float_t         mu2_isGenMatched;
   Float_t         ele1_pt;
   Float_t         ele1_conept;
   Float_t         ele1_eta;
   Float_t         ele1_phi;
   Float_t         ele1_E;
   Int_t           ele1_charge;
   Float_t         ele1_miniRelIso;
   Float_t         ele1_PFRelIso04;
   Int_t           ele1_jetNDauChargedMVASel;
   Float_t         ele1_jetPtRel;
   Float_t         ele1_jetRelIso;
   Float_t         ele1_jetDeepJet;
   Float_t         ele1_sip3D;
   Float_t         ele1_dxy;
   Float_t         ele1_dxyAbs;
   Float_t         ele1_dz;
   Float_t         ele1_ntMVAeleID;
   Float_t         ele1_leptonMVA;
   Int_t           ele1_passesConversionVeto;
   Int_t           ele1_nMissingHits;
   Float_t         ele1_sigmaEtaEta;
   Float_t         ele1_HoE;
   Float_t         ele1_OoEminusOoP;
   Int_t           ele1_isfakeablesel;
   Int_t           ele1_ismvasel;
   Int_t           ele1_isGenMatched;
   Float_t         ele2_pt;
   Float_t         ele2_conept;
   Float_t         ele2_eta;
   Float_t         ele2_phi;
   Float_t         ele2_E;
   Int_t           ele2_charge;
   Float_t         ele2_miniRelIso;
   Float_t         ele2_PFRelIso04;
   Int_t           ele2_jetNDauChargedMVASel;
   Float_t         ele2_jetPtRel;
   Float_t         ele2_jetRelIso;
   Float_t         ele2_jetDeepJet;
   Float_t         ele2_sip3D;
   Float_t         ele2_dxy;
   Float_t         ele2_dxyAbs;
   Float_t         ele2_dz;
   Float_t         ele2_ntMVAeleID;
   Float_t         ele2_leptonMVA;
   Int_t           ele2_passesConversionVeto;
   Int_t           ele2_nMissingHits;
   Float_t         ele2_sigmaEtaEta;
   Float_t         ele2_HoE;
   Float_t         ele2_OoEminusOoP;
   Int_t           ele2_isfakeablesel;
   Int_t           ele2_ismvasel;
   Int_t           ele2_isGenMatched;
   Float_t         ak4Jet1_pt;
   Float_t         ak4Jet1_eta;
   Float_t         ak4Jet1_phi;
   Float_t         ak4Jet1_E;
   Float_t         ak4Jet1_CSV;
   Float_t         ak4Jet1_btagSF;
   Float_t         ak4Jet2_pt;
   Float_t         ak4Jet2_eta;
   Float_t         ak4Jet2_phi;
   Float_t         ak4Jet2_E;
   Float_t         ak4Jet2_CSV;
   Float_t         ak4Jet2_btagSF;
   Float_t         ak4Jet3_pt;
   Float_t         ak4Jet3_eta;
   Float_t         ak4Jet3_phi;
   Float_t         ak4Jet3_E;
   Float_t         ak4Jet3_CSV;
   Float_t         ak4Jet3_btagSF;
   Float_t         ak4Jet4_pt;
   Float_t         ak4Jet4_eta;
   Float_t         ak4Jet4_phi;
   Float_t         ak4Jet4_E;
   Float_t         ak4Jet4_CSV;
   Float_t         ak4Jet4_btagSF;
   Float_t         ak4JetVBF1_pt;
   Float_t         ak4JetVBF1_eta;
   Float_t         ak4JetVBF1_phi;
   Float_t         ak4JetVBF1_E;
   Float_t         ak4JetVBF1_CSV;
   Float_t         ak4JetVBF1_btagSF;
   Float_t         ak4JetVBF2_pt;
   Float_t         ak4JetVBF2_eta;
   Float_t         ak4JetVBF2_phi;
   Float_t         ak4JetVBF2_E;
   Float_t         ak4JetVBF2_CSV;
   Float_t         ak4JetVBF2_btagSF;
   Float_t         ak8Jet1_pt;
   Float_t         ak8Jet1_eta;
   Float_t         ak8Jet1_phi;
   Float_t         ak8Jet1_E;
   Float_t         ak8Jet1_msoftdrop;
   Float_t         ak8Jet1_tau1;
   Float_t         ak8Jet1_tau2;
   Float_t         ak8Jet1_subjet0_pt;
   Float_t         ak8Jet1_subjet0_eta;
   Float_t         ak8Jet1_subjet0_phi;
   Float_t         ak8Jet1_subjet0_CSV;
   Float_t         ak8Jet1_subjet1_pt;
   Float_t         ak8Jet1_subjet1_eta;
   Float_t         ak8Jet1_subjet1_phi;
   Float_t         ak8Jet1_subjet1_CSV;
   Float_t         ak8Jet2_pt;
   Float_t         ak8Jet2_eta;
   Float_t         ak8Jet2_phi;
   Float_t         ak8Jet2_E;
   Float_t         ak8Jet2_msoftdrop;
   Float_t         ak8Jet2_tau1;
   Float_t         ak8Jet2_tau2;
   Float_t         ak8Jet2_subjet0_pt;
   Float_t         ak8Jet2_subjet0_eta;
   Float_t         ak8Jet2_subjet0_phi;
   Float_t         ak8Jet2_subjet0_CSV;
   Float_t         ak8Jet2_subjet1_pt;
   Float_t         ak8Jet2_subjet1_eta;
   Float_t         ak8Jet2_subjet1_phi;
   Float_t         ak8Jet2_subjet1_CSV;
   Float_t         PFMET;
   Float_t         PFMETphi;
   Float_t         HME;
   Float_t         PU_weight;
   Float_t         PU_jetID_SF;
   Float_t         MC_weight;
   Float_t         topPt_wgt;
   Float_t         btag_SF;
   Float_t         trigger_SF;
   Float_t         lepton_IDSF;
   Float_t         lepton_IDSF_recoToLoose;
   Float_t         lepton_IDSF_looseToTight;
   Float_t         L1prefire;
   Float_t         fakeRate;
   Float_t         vbf_m_jj;
   Float_t         vbf_dEta_jj;
   Float_t         trigger_SF_down;
   Float_t         trigger_SF_up;
   Float_t         lepton1_{branchname}_nominal;
   Float_t         lepton1_{branchname}_up;
   Float_t         lepton1_{branchname}_down;
   Float_t         lepton2_{branchname}_nominal;
   Float_t         lepton2_{branchname}_up;
   Float_t         lepton2_{branchname}_down;
   Float_t         genMET_pT;
   Float_t         genMET_phi;
   Float_t         met_covxx;
   Float_t         met_covyy;
   Float_t         met_covxy;
   Int_t           findAllGen;
   Int_t           finalxtohhtobbww;
   Int_t           is_genboosted;
   Int_t           is_genDL;
   Int_t           is_genSL;
   Int_t           is_genee;
   Int_t           is_genem;
   Int_t           is_genmm;
   Float_t         genhtobb_pt;
   Float_t         genhtobb_eta;
   Float_t         genhtobb_phi;
   Float_t         genhtobb_mass;
   Int_t           genhtobb_index;
   Float_t         genhtobb_gendR;
   Int_t           genhtobb_matchidx;
   Int_t           genhtobb_genPartIdxMother;
   Int_t           genhtobb_pdgId;
   Int_t           genhtobb_statusFlags;
   Float_t         genbjet1_pt;
   Float_t         genbjet1_eta;
   Float_t         genbjet1_phi;
   Float_t         genbjet1_mass;
   Int_t           genbjet1_index;
   Float_t         genbjet1_gendR;
   Int_t           genbjet1_matchidx;
   Int_t           genbjet1_hadronFlavour;
   Int_t           genbjet1_partonFlavour;
   Float_t         genbjet1_partondR;
   Float_t         genq2_pt;
   Float_t         genq2_eta;
   Float_t         genq2_phi;
   Float_t         genq2_mass;
   Int_t           genq2_index;
   Float_t         genq2_gendR;
   Int_t           genq2_matchidx;
   Int_t           genq2_genPartIdxMother;
   Int_t           genq2_pdgId;
   Int_t           genq2_statusFlags;
   Float_t         genbjet2nu_pt;
   Float_t         genbjet2nu_eta;
   Float_t         genbjet2nu_phi;
   Float_t         genbjet2nu_mass;
   Int_t           genbjet2nu_index;
   Float_t         genbjet2nu_gendR;
   Int_t           genbjet2nu_matchidx;
   Int_t           genbjet2nu_genPartIdxMother;
   Int_t           genbjet2nu_pdgId;
   Int_t           genbjet2nu_statusFlags;
   Float_t         genak8jet_pt;
   Float_t         genak8jet_eta;
   Float_t         genak8jet_phi;
   Float_t         genak8jet_mass;
   Int_t           genak8jet_index;
   Float_t         genak8jet_gendR;
   Int_t           genak8jet_matchidx;
   Int_t           genak8jet_hadronFlavour;
   Int_t           genak8jet_partonFlavour;
   Float_t         genak8jet_partondR;
   Float_t         genbjet1nu_pt;
   Float_t         genbjet1nu_eta;
   Float_t         genbjet1nu_phi;
   Float_t         genbjet1nu_mass;
   Int_t           genbjet1nu_index;
   Float_t         genbjet1nu_gendR;
   Int_t           genbjet1nu_matchidx;
   Int_t           genbjet1nu_genPartIdxMother;
   Int_t           genbjet1nu_pdgId;
   Int_t           genbjet1nu_statusFlags;
   Float_t         gennu1_pt;
   Float_t         gennu1_eta;
   Float_t         gennu1_phi;
   Float_t         gennu1_mass;
   Int_t           gennu1_index;
   Float_t         gennu1_gendR;
   Int_t           gennu1_matchidx;
   Int_t           gennu1_genPartIdxMother;
   Int_t           gennu1_pdgId;
   Int_t           gennu1_statusFlags;
   Float_t         genbjet2_pt;
   Float_t         genbjet2_eta;
   Float_t         genbjet2_phi;
   Float_t         genbjet2_mass;
   Int_t           genbjet2_index;
   Float_t         genbjet2_gendR;
   Int_t           genbjet2_matchidx;
   Int_t           genbjet2_hadronFlavour;
   Int_t           genbjet2_partonFlavour;
   Float_t         genbjet2_partondR;
   Float_t         genb2_pt;
   Float_t         genb2_eta;
   Float_t         genb2_phi;
   Float_t         genb2_mass;
   Int_t           genb2_index;
   Float_t         genb2_gendR;
   Int_t           genb2_matchidx;
   Int_t           genb2_genPartIdxMother;
   Int_t           genb2_pdgId;
   Int_t           genb2_statusFlags;
   Float_t         gennu2_pt;
   Float_t         gennu2_eta;
   Float_t         gennu2_phi;
   Float_t         gennu2_mass;
   Int_t           gennu2_index;
   Float_t         gennu2_gendR;
   Int_t           gennu2_matchidx;
   Int_t           gennu2_genPartIdxMother;
   Int_t           gennu2_pdgId;
   Int_t           gennu2_statusFlags;
   Float_t         genw2_pt;
   Float_t         genw2_eta;
   Float_t         genw2_phi;
   Float_t         genw2_mass;
   Int_t           genw2_index;
   Float_t         genw2_gendR;
   Int_t           genw2_matchidx;
   Int_t           genw2_genPartIdxMother;
   Int_t           genw2_pdgId;
   Int_t           genw2_statusFlags;
   Float_t         genqjet1_pt;
   Float_t         genqjet1_eta;
   Float_t         genqjet1_phi;
   Float_t         genqjet1_mass;
   Int_t           genqjet1_index;
   Float_t         genqjet1_gendR;
   Int_t           genqjet1_matchidx;
   Int_t           genqjet1_hadronFlavour;
   Int_t           genqjet1_partonFlavour;
   Float_t         genqjet1_partondR;
   Float_t         genqjet2_pt;
   Float_t         genqjet2_eta;
   Float_t         genqjet2_phi;
   Float_t         genqjet2_mass;
   Int_t           genqjet2_index;
   Float_t         genqjet2_gendR;
   Int_t           genqjet2_matchidx;
   Int_t           genqjet2_hadronFlavour;
   Int_t           genqjet2_partonFlavour;
   Float_t         genqjet2_partondR;
   Float_t         genw1_pt;
   Float_t         genw1_eta;
   Float_t         genw1_phi;
   Float_t         genw1_mass;
   Int_t           genw1_index;
   Float_t         genw1_gendR;
   Int_t           genw1_matchidx;
   Int_t           genw1_genPartIdxMother;
   Int_t           genw1_pdgId;
   Int_t           genw1_statusFlags;
   Float_t         genl2_pt;
   Float_t         genl2_eta;
   Float_t         genl2_phi;
   Float_t         genl2_mass;
   Int_t           genl2_index;
   Float_t         genl2_gendR;
   Int_t           genl2_matchidx;
   Int_t           genl2_genPartIdxMother;
   Int_t           genl2_pdgId;
   Int_t           genl2_statusFlags;
   Float_t         genl1_pt;
   Float_t         genl1_eta;
   Float_t         genl1_phi;
   Float_t         genl1_mass;
   Int_t           genl1_index;
   Float_t         genl1_gendR;
   Int_t           genl1_matchidx;
   Int_t           genl1_genPartIdxMother;
   Int_t           genl1_pdgId;
   Int_t           genl1_statusFlags;
   Float_t         genxtohh_pt;
   Float_t         genxtohh_eta;
   Float_t         genxtohh_phi;
   Float_t         genxtohh_mass;
   Int_t           genxtohh_index;
   Float_t         genxtohh_gendR;
   Int_t           genxtohh_matchidx;
   Int_t           genxtohh_genPartIdxMother;
   Int_t           genxtohh_pdgId;
   Int_t           genxtohh_statusFlags;
   Float_t         genhtoww_pt;
   Float_t         genhtoww_eta;
   Float_t         genhtoww_phi;
   Float_t         genhtoww_mass;
   Int_t           genhtoww_index;
   Float_t         genhtoww_gendR;
   Int_t           genhtoww_matchidx;
   Int_t           genhtoww_genPartIdxMother;
   Int_t           genhtoww_pdgId;
   Int_t           genhtoww_statusFlags;
   Float_t         genb1_pt;
   Float_t         genb1_eta;
   Float_t         genb1_phi;
   Float_t         genb1_mass;
   Int_t           genb1_index;
   Float_t         genb1_gendR;
   Int_t           genb1_matchidx;
   Int_t           genb1_genPartIdxMother;
   Int_t           genb1_pdgId;
   Int_t           genb1_statusFlags;
   Float_t         genq1_pt;
   Float_t         genq1_eta;
   Float_t         genq1_phi;
   Float_t         genq1_mass;
   Int_t           genq1_index;
   Float_t         genq1_gendR;
   Int_t           genq1_matchidx;
   Int_t           genq1_genPartIdxMother;
   Int_t           genq1_pdgId;
   Int_t           genq1_statusFlags;

   // List of branches
   TBranch        *b_puWeight;   //!
   TBranch        *b_puWeightUp;   //!
   TBranch        *b_puWeightDown;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_btagSF_deepjet_M_down;   //!
   TBranch        *b_Jet_btagSF_deepjet_M;   //!
   TBranch        *b_Jet_btagSF_deepjet_M_up;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_hf;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_cferr1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_jes;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_cferr2;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_lf;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_lf;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_cferr1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_lfstats1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_lfstats2;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_hfstats1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_hfstats2;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_lfstats2;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_hf;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_lfstats1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_jes;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_hfstats2;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_down_hfstats1;   //!
   TBranch        *b_Jet_btagSF_deepjet_shape_up_cferr2;   //!
   TBranch        *b_nGenHiggs;   //!
   TBranch        *b_GenHiggs_status;   //!
   TBranch        *b_GenHiggs_phi;   //!
   TBranch        *b_GenHiggs_eta;   //!
   TBranch        *b_GenHiggs_mass;   //!
   TBranch        *b_GenHiggs_idx;   //!
   TBranch        *b_GenHiggs_pt;   //!
   TBranch        *b_GenHiggs_statusFlags;   //!
   TBranch        *b_GenHiggs_pdgId;   //!
   TBranch        *b_GenHiggs_charge;   //!
   TBranch        *b_nGenLepFromW1FromHiggs;   //!
   TBranch        *b_GenLepFromW1FromHiggs_status;   //!
   TBranch        *b_GenLepFromW1FromHiggs_phi;   //!
   TBranch        *b_GenLepFromW1FromHiggs_eta;   //!
   TBranch        *b_GenLepFromW1FromHiggs_mass;   //!
   TBranch        *b_GenLepFromW1FromHiggs_idx;   //!
   TBranch        *b_GenLepFromW1FromHiggs_pt;   //!
   TBranch        *b_GenLepFromW1FromHiggs_statusFlags;   //!
   TBranch        *b_GenLepFromW1FromHiggs_pdgId;   //!
   TBranch        *b_GenLepFromW1FromHiggs_charge;   //!
   TBranch        *b_nGenX;   //!
   TBranch        *b_GenX_status;   //!
   TBranch        *b_GenX_phi;   //!
   TBranch        *b_GenX_eta;   //!
   TBranch        *b_GenX_mass;   //!
   TBranch        *b_GenX_idx;   //!
   TBranch        *b_GenX_pt;   //!
   TBranch        *b_GenX_statusFlags;   //!
   TBranch        *b_GenX_pdgId;   //!
   TBranch        *b_GenX_charge;   //!
   TBranch        *b_nGenLepFromW2FromHiggs;   //!
   TBranch        *b_GenLepFromW2FromHiggs_status;   //!
   TBranch        *b_GenLepFromW2FromHiggs_phi;   //!
   TBranch        *b_GenLepFromW2FromHiggs_eta;   //!
   TBranch        *b_GenLepFromW2FromHiggs_mass;   //!
   TBranch        *b_GenLepFromW2FromHiggs_idx;   //!
   TBranch        *b_GenLepFromW2FromHiggs_pt;   //!
   TBranch        *b_GenLepFromW2FromHiggs_statusFlags;   //!
   TBranch        *b_GenLepFromW2FromHiggs_pdgId;   //!
   TBranch        *b_GenLepFromW2FromHiggs_charge;   //!
   TBranch        *b_nGenW1FromHiggs;   //!
   TBranch        *b_GenW1FromHiggs_status;   //!
   TBranch        *b_GenW1FromHiggs_phi;   //!
   TBranch        *b_GenW1FromHiggs_eta;   //!
   TBranch        *b_GenW1FromHiggs_mass;   //!
   TBranch        *b_GenW1FromHiggs_idx;   //!
   TBranch        *b_GenW1FromHiggs_pt;   //!
   TBranch        *b_GenW1FromHiggs_statusFlags;   //!
   TBranch        *b_GenW1FromHiggs_pdgId;   //!
   TBranch        *b_GenW1FromHiggs_charge;   //!
   TBranch        *b_nGenBQuarkFromHiggs;   //!
   TBranch        *b_GenBQuarkFromHiggs_status;   //!
   TBranch        *b_GenBQuarkFromHiggs_phi;   //!
   TBranch        *b_GenBQuarkFromHiggs_eta;   //!
   TBranch        *b_GenBQuarkFromHiggs_mass;   //!
   TBranch        *b_GenBQuarkFromHiggs_idx;   //!
   TBranch        *b_GenBQuarkFromHiggs_pt;   //!
   TBranch        *b_GenBQuarkFromHiggs_statusFlags;   //!
   TBranch        *b_GenBQuarkFromHiggs_pdgId;   //!
   TBranch        *b_GenBQuarkFromHiggs_charge;   //!
   TBranch        *b_nGenW2FromHiggs;   //!
   TBranch        *b_GenW2FromHiggs_status;   //!
   TBranch        *b_GenW2FromHiggs_phi;   //!
   TBranch        *b_GenW2FromHiggs_eta;   //!
   TBranch        *b_GenW2FromHiggs_mass;   //!
   TBranch        *b_GenW2FromHiggs_idx;   //!
   TBranch        *b_GenW2FromHiggs_pt;   //!
   TBranch        *b_GenW2FromHiggs_statusFlags;   //!
   TBranch        *b_GenW2FromHiggs_pdgId;   //!
   TBranch        *b_GenW2FromHiggs_charge;   //!
   TBranch        *b_nGenNuFromW2FromHiggs;   //!
   TBranch        *b_GenNuFromW2FromHiggs_status;   //!
   TBranch        *b_GenNuFromW2FromHiggs_phi;   //!
   TBranch        *b_GenNuFromW2FromHiggs_eta;   //!
   TBranch        *b_GenNuFromW2FromHiggs_mass;   //!
   TBranch        *b_GenNuFromW2FromHiggs_idx;   //!
   TBranch        *b_GenNuFromW2FromHiggs_pt;   //!
   TBranch        *b_GenNuFromW2FromHiggs_statusFlags;   //!
   TBranch        *b_GenNuFromW2FromHiggs_pdgId;   //!
   TBranch        *b_GenNuFromW2FromHiggs_charge;   //!
   TBranch        *b_nGenQuarkFromW2FromHiggs;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_status;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_phi;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_eta;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_mass;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_idx;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_pt;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_statusFlags;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_pdgId;   //!
   TBranch        *b_GenQuarkFromW2FromHiggs_charge;   //!
   TBranch        *b_nGenNuFromW1FromHiggs;   //!
   TBranch        *b_GenNuFromW1FromHiggs_status;   //!
   TBranch        *b_GenNuFromW1FromHiggs_phi;   //!
   TBranch        *b_GenNuFromW1FromHiggs_eta;   //!
   TBranch        *b_GenNuFromW1FromHiggs_mass;   //!
   TBranch        *b_GenNuFromW1FromHiggs_idx;   //!
   TBranch        *b_GenNuFromW1FromHiggs_pt;   //!
   TBranch        *b_GenNuFromW1FromHiggs_statusFlags;   //!
   TBranch        *b_GenNuFromW1FromHiggs_pdgId;   //!
   TBranch        *b_GenNuFromW1FromHiggs_charge;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_run;   //!
   TBranch        *b_n_presel_mu;   //!
   TBranch        *b_n_fakeablesel_mu;   //!
   TBranch        *b_n_mvasel_mu;   //!
   TBranch        *b_n_presel_ele;   //!
   TBranch        *b_n_fakeablesel_ele;   //!
   TBranch        *b_n_mvasel_ele;   //!
   TBranch        *b_n_presel_ak4Jet;   //!
   TBranch        *b_n_presel_ak4JetVBF;   //!
   TBranch        *b_n_presel_ak8Jet;   //!
   TBranch        *b_n_presel_ak8lsJet;   //!
   TBranch        *b_n_loose_ak4BJet;   //!
   TBranch        *b_n_medium_ak4BJet;   //!
   TBranch        *b_is_ee;   //!
   TBranch        *b_is_mm;   //!
   TBranch        *b_is_em;   //!
   TBranch        *b_is_boosted;   //!
   TBranch        *b_is_semiboosted;   //!
   TBranch        *b_is_resolved;   //!
   TBranch        *b_mu1_pt;   //!
   TBranch        *b_mu1_conept;   //!
   TBranch        *b_mu1_eta;   //!
   TBranch        *b_mu1_phi;   //!
   TBranch        *b_mu1_E;   //!
   TBranch        *b_mu1_charge;   //!
   TBranch        *b_mu1_miniRelIso;   //!
   TBranch        *b_mu1_PFRelIso04;   //!
   TBranch        *b_mu1_jetNDauChargedMVASel;   //!
   TBranch        *b_mu1_jetPtRel;   //!
   TBranch        *b_mu1_jetRelIso;   //!
   TBranch        *b_mu1_jetDeepJet;   //!
   TBranch        *b_mu1_sip3D;   //!
   TBranch        *b_mu1_dxy;   //!
   TBranch        *b_mu1_dxyAbs;   //!
   TBranch        *b_mu1_dz;   //!
   TBranch        *b_mu1_segmentCompatibility;   //!
   TBranch        *b_mu1_leptonMVA;   //!
   TBranch        *b_mu1_mediumID;   //!
   TBranch        *b_mu1_dpt_div_pt;   //!
   TBranch        *b_mu1_isfakeablesel;   //!
   TBranch        *b_mu1_ismvasel;   //!
   TBranch        *b_mu1_isGenMatched;   //!
   TBranch        *b_mu2_pt;   //!
   TBranch        *b_mu2_conept;   //!
   TBranch        *b_mu2_eta;   //!
   TBranch        *b_mu2_phi;   //!
   TBranch        *b_mu2_E;   //!
   TBranch        *b_mu2_charge;   //!
   TBranch        *b_mu2_miniRelIso;   //!
   TBranch        *b_mu2_PFRelIso04;   //!
   TBranch        *b_mu2_jetNDauChargedMVASel;   //!
   TBranch        *b_mu2_jetPtRel;   //!
   TBranch        *b_mu2_jetRelIso;   //!
   TBranch        *b_mu2_jetDeepJet;   //!
   TBranch        *b_mu2_sip3D;   //!
   TBranch        *b_mu2_dxy;   //!
   TBranch        *b_mu2_dxyAbs;   //!
   TBranch        *b_mu2_dz;   //!
   TBranch        *b_mu2_segmentCompatibility;   //!
   TBranch        *b_mu2_leptonMVA;   //!
   TBranch        *b_mu2_mediumID;   //!
   TBranch        *b_mu2_dpt_div_pt;   //!
   TBranch        *b_mu2_isfakeablesel;   //!
   TBranch        *b_mu2_ismvasel;   //!
   TBranch        *b_mu2_isGenMatched;   //!
   TBranch        *b_ele1_pt;   //!
   TBranch        *b_ele1_conept;   //!
   TBranch        *b_ele1_eta;   //!
   TBranch        *b_ele1_phi;   //!
   TBranch        *b_ele1_E;   //!
   TBranch        *b_ele1_charge;   //!
   TBranch        *b_ele1_miniRelIso;   //!
   TBranch        *b_ele1_PFRelIso04;   //!
   TBranch        *b_ele1_jetNDauChargedMVASel;   //!
   TBranch        *b_ele1_jetPtRel;   //!
   TBranch        *b_ele1_jetRelIso;   //!
   TBranch        *b_ele1_jetDeepJet;   //!
   TBranch        *b_ele1_sip3D;   //!
   TBranch        *b_ele1_dxy;   //!
   TBranch        *b_ele1_dxyAbs;   //!
   TBranch        *b_ele1_dz;   //!
   TBranch        *b_ele1_ntMVAeleID;   //!
   TBranch        *b_ele1_leptonMVA;   //!
   TBranch        *b_ele1_passesConversionVeto;   //!
   TBranch        *b_ele1_nMissingHits;   //!
   TBranch        *b_ele1_sigmaEtaEta;   //!
   TBranch        *b_ele1_HoE;   //!
   TBranch        *b_ele1_OoEminusOoP;   //!
   TBranch        *b_ele1_isfakeablesel;   //!
   TBranch        *b_ele1_ismvasel;   //!
   TBranch        *b_ele1_isGenMatched;   //!
   TBranch        *b_ele2_pt;   //!
   TBranch        *b_ele2_conept;   //!
   TBranch        *b_ele2_eta;   //!
   TBranch        *b_ele2_phi;   //!
   TBranch        *b_ele2_E;   //!
   TBranch        *b_ele2_charge;   //!
   TBranch        *b_ele2_miniRelIso;   //!
   TBranch        *b_ele2_PFRelIso04;   //!
   TBranch        *b_ele2_jetNDauChargedMVASel;   //!
   TBranch        *b_ele2_jetPtRel;   //!
   TBranch        *b_ele2_jetRelIso;   //!
   TBranch        *b_ele2_jetDeepJet;   //!
   TBranch        *b_ele2_sip3D;   //!
   TBranch        *b_ele2_dxy;   //!
   TBranch        *b_ele2_dxyAbs;   //!
   TBranch        *b_ele2_dz;   //!
   TBranch        *b_ele2_ntMVAeleID;   //!
   TBranch        *b_ele2_leptonMVA;   //!
   TBranch        *b_ele2_passesConversionVeto;   //!
   TBranch        *b_ele2_nMissingHits;   //!
   TBranch        *b_ele2_sigmaEtaEta;   //!
   TBranch        *b_ele2_HoE;   //!
   TBranch        *b_ele2_OoEminusOoP;   //!
   TBranch        *b_ele2_isfakeablesel;   //!
   TBranch        *b_ele2_ismvasel;   //!
   TBranch        *b_ele2_isGenMatched;   //!
   TBranch        *b_ak4Jet1_pt;   //!
   TBranch        *b_ak4Jet1_eta;   //!
   TBranch        *b_ak4Jet1_phi;   //!
   TBranch        *b_ak4Jet1_E;   //!
   TBranch        *b_ak4Jet1_CSV;   //!
   TBranch        *b_ak4Jet1_btagSF;   //!
   TBranch        *b_ak4Jet2_pt;   //!
   TBranch        *b_ak4Jet2_eta;   //!
   TBranch        *b_ak4Jet2_phi;   //!
   TBranch        *b_ak4Jet2_E;   //!
   TBranch        *b_ak4Jet2_CSV;   //!
   TBranch        *b_ak4Jet2_btagSF;   //!
   TBranch        *b_ak4Jet3_pt;   //!
   TBranch        *b_ak4Jet3_eta;   //!
   TBranch        *b_ak4Jet3_phi;   //!
   TBranch        *b_ak4Jet3_E;   //!
   TBranch        *b_ak4Jet3_CSV;   //!
   TBranch        *b_ak4Jet3_btagSF;   //!
   TBranch        *b_ak4Jet4_pt;   //!
   TBranch        *b_ak4Jet4_eta;   //!
   TBranch        *b_ak4Jet4_phi;   //!
   TBranch        *b_ak4Jet4_E;   //!
   TBranch        *b_ak4Jet4_CSV;   //!
   TBranch        *b_ak4Jet4_btagSF;   //!
   TBranch        *b_ak4JetVBF1_pt;   //!
   TBranch        *b_ak4JetVBF1_eta;   //!
   TBranch        *b_ak4JetVBF1_phi;   //!
   TBranch        *b_ak4JetVBF1_E;   //!
   TBranch        *b_ak4JetVBF1_CSV;   //!
   TBranch        *b_ak4JetVBF1_btagSF;   //!
   TBranch        *b_ak4JetVBF2_pt;   //!
   TBranch        *b_ak4JetVBF2_eta;   //!
   TBranch        *b_ak4JetVBF2_phi;   //!
   TBranch        *b_ak4JetVBF2_E;   //!
   TBranch        *b_ak4JetVBF2_CSV;   //!
   TBranch        *b_ak4JetVBF2_btagSF;   //!
   TBranch        *b_ak8Jet1_pt;   //!
   TBranch        *b_ak8Jet1_eta;   //!
   TBranch        *b_ak8Jet1_phi;   //!
   TBranch        *b_ak8Jet1_E;   //!
   TBranch        *b_ak8Jet1_msoftdrop;   //!
   TBranch        *b_ak8Jet1_tau1;   //!
   TBranch        *b_ak8Jet1_tau2;   //!
   TBranch        *b_ak8Jet1_subjet0_pt;   //!
   TBranch        *b_ak8Jet1_subjet0_eta;   //!
   TBranch        *b_ak8Jet1_subjet0_phi;   //!
   TBranch        *b_ak8Jet1_subjet0_CSV;   //!
   TBranch        *b_ak8Jet1_subjet1_pt;   //!
   TBranch        *b_ak8Jet1_subjet1_eta;   //!
   TBranch        *b_ak8Jet1_subjet1_phi;   //!
   TBranch        *b_ak8Jet1_subjet1_CSV;   //!
   TBranch        *b_ak8Jet2_pt;   //!
   TBranch        *b_ak8Jet2_eta;   //!
   TBranch        *b_ak8Jet2_phi;   //!
   TBranch        *b_ak8Jet2_E;   //!
   TBranch        *b_ak8Jet2_msoftdrop;   //!
   TBranch        *b_ak8Jet2_tau1;   //!
   TBranch        *b_ak8Jet2_tau2;   //!
   TBranch        *b_ak8Jet2_subjet0_pt;   //!
   TBranch        *b_ak8Jet2_subjet0_eta;   //!
   TBranch        *b_ak8Jet2_subjet0_phi;   //!
   TBranch        *b_ak8Jet2_subjet0_CSV;   //!
   TBranch        *b_ak8Jet2_subjet1_pt;   //!
   TBranch        *b_ak8Jet2_subjet1_eta;   //!
   TBranch        *b_ak8Jet2_subjet1_phi;   //!
   TBranch        *b_ak8Jet2_subjet1_CSV;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_PFMETphi;   //!
   TBranch        *b_HME;   //!
   TBranch        *b_PU_weight;   //!
   TBranch        *b_PU_jetID_SF;   //!
   TBranch        *b_MC_weight;   //!
   TBranch        *b_topPt_wgt;   //!
   TBranch        *b_btag_SF;   //!
   TBranch        *b_trigger_SF;   //!
   TBranch        *b_lepton_IDSF;   //!
   TBranch        *b_lepton_IDSF_recoToLoose;   //!
   TBranch        *b_lepton_IDSF_looseToTight;   //!
   TBranch        *b_L1prefire;   //!
   TBranch        *b_fakeRate;   //!
   TBranch        *b_vbf_m_jj;   //!
   TBranch        *b_vbf_dEta_jj;   //!
   TBranch        *b_trigger_SF_down;   //!
   TBranch        *b_trigger_SF_up;   //!
   TBranch        *b_lepton1_{branchname}_nominal;   //!
   TBranch        *b_lepton1_{branchname}_up;   //!
   TBranch        *b_lepton1_{branchname}_down;   //!
   TBranch        *b_lepton2_{branchname}_nominal;   //!
   TBranch        *b_lepton2_{branchname}_up;   //!
   TBranch        *b_lepton2_{branchname}_down;   //!
   TBranch        *b_genMET_pT;   //!
   TBranch        *b_genMET_phi;   //!
   TBranch        *b_met_covxx;   //!
   TBranch        *b_met_covyy;   //!
   TBranch        *b_met_covxy;   //!
   TBranch        *b_findAllGen;   //!
   TBranch        *b_finalxtohhtobbww;   //!
   TBranch        *b_is_genboosted;   //!
   TBranch        *b_is_genDL;   //!
   TBranch        *b_is_genSL;   //!
   TBranch        *b_is_genee;   //!
   TBranch        *b_is_genem;   //!
   TBranch        *b_is_genmm;   //!
   TBranch        *b_genhtobb_pt;   //!
   TBranch        *b_genhtobb_eta;   //!
   TBranch        *b_genhtobb_phi;   //!
   TBranch        *b_genhtobb_mass;   //!
   TBranch        *b_genhtobb_index;   //!
   TBranch        *b_genhtobb_gendR;   //!
   TBranch        *b_genhtobb_matchidx;   //!
   TBranch        *b_genhtobb_genPartIdxMother;   //!
   TBranch        *b_genhtobb_pdgId;   //!
   TBranch        *b_genhtobb_statusFlags;   //!
   TBranch        *b_genbjet1_pt;   //!
   TBranch        *b_genbjet1_eta;   //!
   TBranch        *b_genbjet1_phi;   //!
   TBranch        *b_genbjet1_mass;   //!
   TBranch        *b_genbjet1_index;   //!
   TBranch        *b_genbjet1_gendR;   //!
   TBranch        *b_genbjet1_matchidx;   //!
   TBranch        *b_genbjet1_hadronFlavour;   //!
   TBranch        *b_genbjet1_partonFlavour;   //!
   TBranch        *b_genbjet1_partondR;   //!
   TBranch        *b_genq2_pt;   //!
   TBranch        *b_genq2_eta;   //!
   TBranch        *b_genq2_phi;   //!
   TBranch        *b_genq2_mass;   //!
   TBranch        *b_genq2_index;   //!
   TBranch        *b_genq2_gendR;   //!
   TBranch        *b_genq2_matchidx;   //!
   TBranch        *b_genq2_genPartIdxMother;   //!
   TBranch        *b_genq2_pdgId;   //!
   TBranch        *b_genq2_statusFlags;   //!
   TBranch        *b_genbjet2nu_pt;   //!
   TBranch        *b_genbjet2nu_eta;   //!
   TBranch        *b_genbjet2nu_phi;   //!
   TBranch        *b_genbjet2nu_mass;   //!
   TBranch        *b_genbjet2nu_index;   //!
   TBranch        *b_genbjet2nu_gendR;   //!
   TBranch        *b_genbjet2nu_matchidx;   //!
   TBranch        *b_genbjet2nu_genPartIdxMother;   //!
   TBranch        *b_genbjet2nu_pdgId;   //!
   TBranch        *b_genbjet2nu_statusFlags;   //!
   TBranch        *b_genak8jet_pt;   //!
   TBranch        *b_genak8jet_eta;   //!
   TBranch        *b_genak8jet_phi;   //!
   TBranch        *b_genak8jet_mass;   //!
   TBranch        *b_genak8jet_index;   //!
   TBranch        *b_genak8jet_gendR;   //!
   TBranch        *b_genak8jet_matchidx;   //!
   TBranch        *b_genak8jet_hadronFlavour;   //!
   TBranch        *b_genak8jet_partonFlavour;   //!
   TBranch        *b_genak8jet_partondR;   //!
   TBranch        *b_genbjet1nu_pt;   //!
   TBranch        *b_genbjet1nu_eta;   //!
   TBranch        *b_genbjet1nu_phi;   //!
   TBranch        *b_genbjet1nu_mass;   //!
   TBranch        *b_genbjet1nu_index;   //!
   TBranch        *b_genbjet1nu_gendR;   //!
   TBranch        *b_genbjet1nu_matchidx;   //!
   TBranch        *b_genbjet1nu_genPartIdxMother;   //!
   TBranch        *b_genbjet1nu_pdgId;   //!
   TBranch        *b_genbjet1nu_statusFlags;   //!
   TBranch        *b_gennu1_pt;   //!
   TBranch        *b_gennu1_eta;   //!
   TBranch        *b_gennu1_phi;   //!
   TBranch        *b_gennu1_mass;   //!
   TBranch        *b_gennu1_index;   //!
   TBranch        *b_gennu1_gendR;   //!
   TBranch        *b_gennu1_matchidx;   //!
   TBranch        *b_gennu1_genPartIdxMother;   //!
   TBranch        *b_gennu1_pdgId;   //!
   TBranch        *b_gennu1_statusFlags;   //!
   TBranch        *b_genbjet2_pt;   //!
   TBranch        *b_genbjet2_eta;   //!
   TBranch        *b_genbjet2_phi;   //!
   TBranch        *b_genbjet2_mass;   //!
   TBranch        *b_genbjet2_index;   //!
   TBranch        *b_genbjet2_gendR;   //!
   TBranch        *b_genbjet2_matchidx;   //!
   TBranch        *b_genbjet2_hadronFlavour;   //!
   TBranch        *b_genbjet2_partonFlavour;   //!
   TBranch        *b_genbjet2_partondR;   //!
   TBranch        *b_genb2_pt;   //!
   TBranch        *b_genb2_eta;   //!
   TBranch        *b_genb2_phi;   //!
   TBranch        *b_genb2_mass;   //!
   TBranch        *b_genb2_index;   //!
   TBranch        *b_genb2_gendR;   //!
   TBranch        *b_genb2_matchidx;   //!
   TBranch        *b_genb2_genPartIdxMother;   //!
   TBranch        *b_genb2_pdgId;   //!
   TBranch        *b_genb2_statusFlags;   //!
   TBranch        *b_gennu2_pt;   //!
   TBranch        *b_gennu2_eta;   //!
   TBranch        *b_gennu2_phi;   //!
   TBranch        *b_gennu2_mass;   //!
   TBranch        *b_gennu2_index;   //!
   TBranch        *b_gennu2_gendR;   //!
   TBranch        *b_gennu2_matchidx;   //!
   TBranch        *b_gennu2_genPartIdxMother;   //!
   TBranch        *b_gennu2_pdgId;   //!
   TBranch        *b_gennu2_statusFlags;   //!
   TBranch        *b_genw2_pt;   //!
   TBranch        *b_genw2_eta;   //!
   TBranch        *b_genw2_phi;   //!
   TBranch        *b_genw2_mass;   //!
   TBranch        *b_genw2_index;   //!
   TBranch        *b_genw2_gendR;   //!
   TBranch        *b_genw2_matchidx;   //!
   TBranch        *b_genw2_genPartIdxMother;   //!
   TBranch        *b_genw2_pdgId;   //!
   TBranch        *b_genw2_statusFlags;   //!
   TBranch        *b_genqjet1_pt;   //!
   TBranch        *b_genqjet1_eta;   //!
   TBranch        *b_genqjet1_phi;   //!
   TBranch        *b_genqjet1_mass;   //!
   TBranch        *b_genqjet1_index;   //!
   TBranch        *b_genqjet1_gendR;   //!
   TBranch        *b_genqjet1_matchidx;   //!
   TBranch        *b_genqjet1_hadronFlavour;   //!
   TBranch        *b_genqjet1_partonFlavour;   //!
   TBranch        *b_genqjet1_partondR;   //!
   TBranch        *b_genqjet2_pt;   //!
   TBranch        *b_genqjet2_eta;   //!
   TBranch        *b_genqjet2_phi;   //!
   TBranch        *b_genqjet2_mass;   //!
   TBranch        *b_genqjet2_index;   //!
   TBranch        *b_genqjet2_gendR;   //!
   TBranch        *b_genqjet2_matchidx;   //!
   TBranch        *b_genqjet2_hadronFlavour;   //!
   TBranch        *b_genqjet2_partonFlavour;   //!
   TBranch        *b_genqjet2_partondR;   //!
   TBranch        *b_genw1_pt;   //!
   TBranch        *b_genw1_eta;   //!
   TBranch        *b_genw1_phi;   //!
   TBranch        *b_genw1_mass;   //!
   TBranch        *b_genw1_index;   //!
   TBranch        *b_genw1_gendR;   //!
   TBranch        *b_genw1_matchidx;   //!
   TBranch        *b_genw1_genPartIdxMother;   //!
   TBranch        *b_genw1_pdgId;   //!
   TBranch        *b_genw1_statusFlags;   //!
   TBranch        *b_genl2_pt;   //!
   TBranch        *b_genl2_eta;   //!
   TBranch        *b_genl2_phi;   //!
   TBranch        *b_genl2_mass;   //!
   TBranch        *b_genl2_index;   //!
   TBranch        *b_genl2_gendR;   //!
   TBranch        *b_genl2_matchidx;   //!
   TBranch        *b_genl2_genPartIdxMother;   //!
   TBranch        *b_genl2_pdgId;   //!
   TBranch        *b_genl2_statusFlags;   //!
   TBranch        *b_genl1_pt;   //!
   TBranch        *b_genl1_eta;   //!
   TBranch        *b_genl1_phi;   //!
   TBranch        *b_genl1_mass;   //!
   TBranch        *b_genl1_index;   //!
   TBranch        *b_genl1_gendR;   //!
   TBranch        *b_genl1_matchidx;   //!
   TBranch        *b_genl1_genPartIdxMother;   //!
   TBranch        *b_genl1_pdgId;   //!
   TBranch        *b_genl1_statusFlags;   //!
   TBranch        *b_genxtohh_pt;   //!
   TBranch        *b_genxtohh_eta;   //!
   TBranch        *b_genxtohh_phi;   //!
   TBranch        *b_genxtohh_mass;   //!
   TBranch        *b_genxtohh_index;   //!
   TBranch        *b_genxtohh_gendR;   //!
   TBranch        *b_genxtohh_matchidx;   //!
   TBranch        *b_genxtohh_genPartIdxMother;   //!
   TBranch        *b_genxtohh_pdgId;   //!
   TBranch        *b_genxtohh_statusFlags;   //!
   TBranch        *b_genhtoww_pt;   //!
   TBranch        *b_genhtoww_eta;   //!
   TBranch        *b_genhtoww_phi;   //!
   TBranch        *b_genhtoww_mass;   //!
   TBranch        *b_genhtoww_index;   //!
   TBranch        *b_genhtoww_gendR;   //!
   TBranch        *b_genhtoww_matchidx;   //!
   TBranch        *b_genhtoww_genPartIdxMother;   //!
   TBranch        *b_genhtoww_pdgId;   //!
   TBranch        *b_genhtoww_statusFlags;   //!
   TBranch        *b_genb1_pt;   //!
   TBranch        *b_genb1_eta;   //!
   TBranch        *b_genb1_phi;   //!
   TBranch        *b_genb1_mass;   //!
   TBranch        *b_genb1_index;   //!
   TBranch        *b_genb1_gendR;   //!
   TBranch        *b_genb1_matchidx;   //!
   TBranch        *b_genb1_genPartIdxMother;   //!
   TBranch        *b_genb1_pdgId;   //!
   TBranch        *b_genb1_statusFlags;   //!
   TBranch        *b_genq1_pt;   //!
   TBranch        *b_genq1_eta;   //!
   TBranch        *b_genq1_phi;   //!
   TBranch        *b_genq1_mass;   //!
   TBranch        *b_genq1_index;   //!
   TBranch        *b_genq1_gendR;   //!
   TBranch        *b_genq1_matchidx;   //!
   TBranch        *b_genq1_genPartIdxMother;   //!
   TBranch        *b_genq1_pdgId;   //!
   TBranch        *b_genq1_statusFlags;   //!

   syncTree(TTree *tree=0);
   virtual ~syncTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef syncTree_cxx
syncTree::syncTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
      }
      f->GetObject("syncTree",tree);

   }
   Init(tree);
}

syncTree::~syncTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t syncTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t syncTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void syncTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   fChain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M_down", Jet_btagSF_deepjet_M_down, &b_Jet_btagSF_deepjet_M_down);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M", Jet_btagSF_deepjet_M, &b_Jet_btagSF_deepjet_M);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_M_up", Jet_btagSF_deepjet_M_up, &b_Jet_btagSF_deepjet_M_up);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_hf", Jet_btagSF_deepjet_shape_down_hf, &b_Jet_btagSF_deepjet_shape_down_hf);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape", Jet_btagSF_deepjet_shape, &b_Jet_btagSF_deepjet_shape);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_cferr1", Jet_btagSF_deepjet_shape_up_cferr1, &b_Jet_btagSF_deepjet_shape_up_cferr1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_jes", Jet_btagSF_deepjet_shape_up_jes, &b_Jet_btagSF_deepjet_shape_up_jes);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_cferr2", Jet_btagSF_deepjet_shape_down_cferr2, &b_Jet_btagSF_deepjet_shape_down_cferr2);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_lf", Jet_btagSF_deepjet_shape_up_lf, &b_Jet_btagSF_deepjet_shape_up_lf);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_lf", Jet_btagSF_deepjet_shape_down_lf, &b_Jet_btagSF_deepjet_shape_down_lf);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_cferr1", Jet_btagSF_deepjet_shape_down_cferr1, &b_Jet_btagSF_deepjet_shape_down_cferr1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_lfstats1", Jet_btagSF_deepjet_shape_up_lfstats1, &b_Jet_btagSF_deepjet_shape_up_lfstats1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_lfstats2", Jet_btagSF_deepjet_shape_up_lfstats2, &b_Jet_btagSF_deepjet_shape_up_lfstats2);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_hfstats1", Jet_btagSF_deepjet_shape_up_hfstats1, &b_Jet_btagSF_deepjet_shape_up_hfstats1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_hfstats2", Jet_btagSF_deepjet_shape_up_hfstats2, &b_Jet_btagSF_deepjet_shape_up_hfstats2);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_lfstats2", Jet_btagSF_deepjet_shape_down_lfstats2, &b_Jet_btagSF_deepjet_shape_down_lfstats2);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_hf", Jet_btagSF_deepjet_shape_up_hf, &b_Jet_btagSF_deepjet_shape_up_hf);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_lfstats1", Jet_btagSF_deepjet_shape_down_lfstats1, &b_Jet_btagSF_deepjet_shape_down_lfstats1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_jes", Jet_btagSF_deepjet_shape_down_jes, &b_Jet_btagSF_deepjet_shape_down_jes);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_hfstats2", Jet_btagSF_deepjet_shape_down_hfstats2, &b_Jet_btagSF_deepjet_shape_down_hfstats2);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_down_hfstats1", Jet_btagSF_deepjet_shape_down_hfstats1, &b_Jet_btagSF_deepjet_shape_down_hfstats1);
   fChain->SetBranchAddress("Jet_btagSF_deepjet_shape_up_cferr2", Jet_btagSF_deepjet_shape_up_cferr2, &b_Jet_btagSF_deepjet_shape_up_cferr2);
   fChain->SetBranchAddress("nGenHiggs", &nGenHiggs, &b_nGenHiggs);
   fChain->SetBranchAddress("GenHiggs_status", GenHiggs_status, &b_GenHiggs_status);
   fChain->SetBranchAddress("GenHiggs_phi", GenHiggs_phi, &b_GenHiggs_phi);
   fChain->SetBranchAddress("GenHiggs_eta", GenHiggs_eta, &b_GenHiggs_eta);
   fChain->SetBranchAddress("GenHiggs_mass", GenHiggs_mass, &b_GenHiggs_mass);
   fChain->SetBranchAddress("GenHiggs_idx", GenHiggs_idx, &b_GenHiggs_idx);
   fChain->SetBranchAddress("GenHiggs_pt", GenHiggs_pt, &b_GenHiggs_pt);
   fChain->SetBranchAddress("GenHiggs_statusFlags", GenHiggs_statusFlags, &b_GenHiggs_statusFlags);
   fChain->SetBranchAddress("GenHiggs_pdgId", GenHiggs_pdgId, &b_GenHiggs_pdgId);
   fChain->SetBranchAddress("GenHiggs_charge", GenHiggs_charge, &b_GenHiggs_charge);
   fChain->SetBranchAddress("nGenLepFromW1FromHiggs", &nGenLepFromW1FromHiggs, &b_nGenLepFromW1FromHiggs);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_status", GenLepFromW1FromHiggs_status, &b_GenLepFromW1FromHiggs_status);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_phi", GenLepFromW1FromHiggs_phi, &b_GenLepFromW1FromHiggs_phi);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_eta", GenLepFromW1FromHiggs_eta, &b_GenLepFromW1FromHiggs_eta);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_mass", GenLepFromW1FromHiggs_mass, &b_GenLepFromW1FromHiggs_mass);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_idx", GenLepFromW1FromHiggs_idx, &b_GenLepFromW1FromHiggs_idx);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_pt", GenLepFromW1FromHiggs_pt, &b_GenLepFromW1FromHiggs_pt);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_statusFlags", GenLepFromW1FromHiggs_statusFlags, &b_GenLepFromW1FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_pdgId", GenLepFromW1FromHiggs_pdgId, &b_GenLepFromW1FromHiggs_pdgId);
   fChain->SetBranchAddress("GenLepFromW1FromHiggs_charge", GenLepFromW1FromHiggs_charge, &b_GenLepFromW1FromHiggs_charge);
   fChain->SetBranchAddress("nGenX", &nGenX, &b_nGenX);
   fChain->SetBranchAddress("GenX_status", GenX_status, &b_GenX_status);
   fChain->SetBranchAddress("GenX_phi", GenX_phi, &b_GenX_phi);
   fChain->SetBranchAddress("GenX_eta", GenX_eta, &b_GenX_eta);
   fChain->SetBranchAddress("GenX_mass", GenX_mass, &b_GenX_mass);
   fChain->SetBranchAddress("GenX_idx", GenX_idx, &b_GenX_idx);
   fChain->SetBranchAddress("GenX_pt", GenX_pt, &b_GenX_pt);
   fChain->SetBranchAddress("GenX_statusFlags", GenX_statusFlags, &b_GenX_statusFlags);
   fChain->SetBranchAddress("GenX_pdgId", GenX_pdgId, &b_GenX_pdgId);
   fChain->SetBranchAddress("GenX_charge", GenX_charge, &b_GenX_charge);
   fChain->SetBranchAddress("nGenLepFromW2FromHiggs", &nGenLepFromW2FromHiggs, &b_nGenLepFromW2FromHiggs);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_status", &GenLepFromW2FromHiggs_status, &b_GenLepFromW2FromHiggs_status);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_phi", &GenLepFromW2FromHiggs_phi, &b_GenLepFromW2FromHiggs_phi);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_eta", &GenLepFromW2FromHiggs_eta, &b_GenLepFromW2FromHiggs_eta);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_mass", &GenLepFromW2FromHiggs_mass, &b_GenLepFromW2FromHiggs_mass);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_idx", &GenLepFromW2FromHiggs_idx, &b_GenLepFromW2FromHiggs_idx);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_pt", &GenLepFromW2FromHiggs_pt, &b_GenLepFromW2FromHiggs_pt);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_statusFlags", &GenLepFromW2FromHiggs_statusFlags, &b_GenLepFromW2FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_pdgId", &GenLepFromW2FromHiggs_pdgId, &b_GenLepFromW2FromHiggs_pdgId);
   fChain->SetBranchAddress("GenLepFromW2FromHiggs_charge", &GenLepFromW2FromHiggs_charge, &b_GenLepFromW2FromHiggs_charge);
   fChain->SetBranchAddress("nGenW1FromHiggs", &nGenW1FromHiggs, &b_nGenW1FromHiggs);
   fChain->SetBranchAddress("GenW1FromHiggs_status", GenW1FromHiggs_status, &b_GenW1FromHiggs_status);
   fChain->SetBranchAddress("GenW1FromHiggs_phi", GenW1FromHiggs_phi, &b_GenW1FromHiggs_phi);
   fChain->SetBranchAddress("GenW1FromHiggs_eta", GenW1FromHiggs_eta, &b_GenW1FromHiggs_eta);
   fChain->SetBranchAddress("GenW1FromHiggs_mass", GenW1FromHiggs_mass, &b_GenW1FromHiggs_mass);
   fChain->SetBranchAddress("GenW1FromHiggs_idx", GenW1FromHiggs_idx, &b_GenW1FromHiggs_idx);
   fChain->SetBranchAddress("GenW1FromHiggs_pt", GenW1FromHiggs_pt, &b_GenW1FromHiggs_pt);
   fChain->SetBranchAddress("GenW1FromHiggs_statusFlags", GenW1FromHiggs_statusFlags, &b_GenW1FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenW1FromHiggs_pdgId", GenW1FromHiggs_pdgId, &b_GenW1FromHiggs_pdgId);
   fChain->SetBranchAddress("GenW1FromHiggs_charge", GenW1FromHiggs_charge, &b_GenW1FromHiggs_charge);
   fChain->SetBranchAddress("nGenBQuarkFromHiggs", &nGenBQuarkFromHiggs, &b_nGenBQuarkFromHiggs);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_status", GenBQuarkFromHiggs_status, &b_GenBQuarkFromHiggs_status);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_phi", GenBQuarkFromHiggs_phi, &b_GenBQuarkFromHiggs_phi);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_eta", GenBQuarkFromHiggs_eta, &b_GenBQuarkFromHiggs_eta);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_mass", GenBQuarkFromHiggs_mass, &b_GenBQuarkFromHiggs_mass);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_idx", GenBQuarkFromHiggs_idx, &b_GenBQuarkFromHiggs_idx);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_pt", GenBQuarkFromHiggs_pt, &b_GenBQuarkFromHiggs_pt);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_statusFlags", GenBQuarkFromHiggs_statusFlags, &b_GenBQuarkFromHiggs_statusFlags);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_pdgId", GenBQuarkFromHiggs_pdgId, &b_GenBQuarkFromHiggs_pdgId);
   fChain->SetBranchAddress("GenBQuarkFromHiggs_charge", GenBQuarkFromHiggs_charge, &b_GenBQuarkFromHiggs_charge);
   fChain->SetBranchAddress("nGenW2FromHiggs", &nGenW2FromHiggs, &b_nGenW2FromHiggs);
   fChain->SetBranchAddress("GenW2FromHiggs_status", GenW2FromHiggs_status, &b_GenW2FromHiggs_status);
   fChain->SetBranchAddress("GenW2FromHiggs_phi", GenW2FromHiggs_phi, &b_GenW2FromHiggs_phi);
   fChain->SetBranchAddress("GenW2FromHiggs_eta", GenW2FromHiggs_eta, &b_GenW2FromHiggs_eta);
   fChain->SetBranchAddress("GenW2FromHiggs_mass", GenW2FromHiggs_mass, &b_GenW2FromHiggs_mass);
   fChain->SetBranchAddress("GenW2FromHiggs_idx", GenW2FromHiggs_idx, &b_GenW2FromHiggs_idx);
   fChain->SetBranchAddress("GenW2FromHiggs_pt", GenW2FromHiggs_pt, &b_GenW2FromHiggs_pt);
   fChain->SetBranchAddress("GenW2FromHiggs_statusFlags", GenW2FromHiggs_statusFlags, &b_GenW2FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenW2FromHiggs_pdgId", GenW2FromHiggs_pdgId, &b_GenW2FromHiggs_pdgId);
   fChain->SetBranchAddress("GenW2FromHiggs_charge", GenW2FromHiggs_charge, &b_GenW2FromHiggs_charge);
   fChain->SetBranchAddress("nGenNuFromW2FromHiggs", &nGenNuFromW2FromHiggs, &b_nGenNuFromW2FromHiggs);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_status", &GenNuFromW2FromHiggs_status, &b_GenNuFromW2FromHiggs_status);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_phi", &GenNuFromW2FromHiggs_phi, &b_GenNuFromW2FromHiggs_phi);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_eta", &GenNuFromW2FromHiggs_eta, &b_GenNuFromW2FromHiggs_eta);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_mass", &GenNuFromW2FromHiggs_mass, &b_GenNuFromW2FromHiggs_mass);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_idx", &GenNuFromW2FromHiggs_idx, &b_GenNuFromW2FromHiggs_idx);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_pt", &GenNuFromW2FromHiggs_pt, &b_GenNuFromW2FromHiggs_pt);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_statusFlags", &GenNuFromW2FromHiggs_statusFlags, &b_GenNuFromW2FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_pdgId", &GenNuFromW2FromHiggs_pdgId, &b_GenNuFromW2FromHiggs_pdgId);
   fChain->SetBranchAddress("GenNuFromW2FromHiggs_charge", &GenNuFromW2FromHiggs_charge, &b_GenNuFromW2FromHiggs_charge);
   fChain->SetBranchAddress("nGenQuarkFromW2FromHiggs", &nGenQuarkFromW2FromHiggs, &b_nGenQuarkFromW2FromHiggs);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_status", GenQuarkFromW2FromHiggs_status, &b_GenQuarkFromW2FromHiggs_status);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_phi", GenQuarkFromW2FromHiggs_phi, &b_GenQuarkFromW2FromHiggs_phi);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_eta", GenQuarkFromW2FromHiggs_eta, &b_GenQuarkFromW2FromHiggs_eta);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_mass", GenQuarkFromW2FromHiggs_mass, &b_GenQuarkFromW2FromHiggs_mass);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_idx", GenQuarkFromW2FromHiggs_idx, &b_GenQuarkFromW2FromHiggs_idx);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_pt", GenQuarkFromW2FromHiggs_pt, &b_GenQuarkFromW2FromHiggs_pt);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_statusFlags", GenQuarkFromW2FromHiggs_statusFlags, &b_GenQuarkFromW2FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_pdgId", GenQuarkFromW2FromHiggs_pdgId, &b_GenQuarkFromW2FromHiggs_pdgId);
   fChain->SetBranchAddress("GenQuarkFromW2FromHiggs_charge", GenQuarkFromW2FromHiggs_charge, &b_GenQuarkFromW2FromHiggs_charge);
   fChain->SetBranchAddress("nGenNuFromW1FromHiggs", &nGenNuFromW1FromHiggs, &b_nGenNuFromW1FromHiggs);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_status", GenNuFromW1FromHiggs_status, &b_GenNuFromW1FromHiggs_status);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_phi", GenNuFromW1FromHiggs_phi, &b_GenNuFromW1FromHiggs_phi);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_eta", GenNuFromW1FromHiggs_eta, &b_GenNuFromW1FromHiggs_eta);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_mass", GenNuFromW1FromHiggs_mass, &b_GenNuFromW1FromHiggs_mass);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_idx", GenNuFromW1FromHiggs_idx, &b_GenNuFromW1FromHiggs_idx);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_pt", GenNuFromW1FromHiggs_pt, &b_GenNuFromW1FromHiggs_pt);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_statusFlags", GenNuFromW1FromHiggs_statusFlags, &b_GenNuFromW1FromHiggs_statusFlags);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_pdgId", GenNuFromW1FromHiggs_pdgId, &b_GenNuFromW1FromHiggs_pdgId);
   fChain->SetBranchAddress("GenNuFromW1FromHiggs_charge", GenNuFromW1FromHiggs_charge, &b_GenNuFromW1FromHiggs_charge);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("n_presel_mu", &n_presel_mu, &b_n_presel_mu);
   fChain->SetBranchAddress("n_fakeablesel_mu", &n_fakeablesel_mu, &b_n_fakeablesel_mu);
   fChain->SetBranchAddress("n_mvasel_mu", &n_mvasel_mu, &b_n_mvasel_mu);
   fChain->SetBranchAddress("n_presel_ele", &n_presel_ele, &b_n_presel_ele);
   fChain->SetBranchAddress("n_fakeablesel_ele", &n_fakeablesel_ele, &b_n_fakeablesel_ele);
   fChain->SetBranchAddress("n_mvasel_ele", &n_mvasel_ele, &b_n_mvasel_ele);
   fChain->SetBranchAddress("n_presel_ak4Jet", &n_presel_ak4Jet, &b_n_presel_ak4Jet);
   fChain->SetBranchAddress("n_presel_ak4JetVBF", &n_presel_ak4JetVBF, &b_n_presel_ak4JetVBF);
   fChain->SetBranchAddress("n_presel_ak8Jet", &n_presel_ak8Jet, &b_n_presel_ak8Jet);
   fChain->SetBranchAddress("n_presel_ak8lsJet", &n_presel_ak8lsJet, &b_n_presel_ak8lsJet);
   fChain->SetBranchAddress("n_loose_ak4BJet", &n_loose_ak4BJet, &b_n_loose_ak4BJet);
   fChain->SetBranchAddress("n_medium_ak4BJet", &n_medium_ak4BJet, &b_n_medium_ak4BJet);
   fChain->SetBranchAddress("is_ee", &is_ee, &b_is_ee);
   fChain->SetBranchAddress("is_mm", &is_mm, &b_is_mm);
   fChain->SetBranchAddress("is_em", &is_em, &b_is_em);
   fChain->SetBranchAddress("is_boosted", &is_boosted, &b_is_boosted);
   fChain->SetBranchAddress("is_semiboosted", &is_semiboosted, &b_is_semiboosted);
   fChain->SetBranchAddress("is_resolved", &is_resolved, &b_is_resolved);
   fChain->SetBranchAddress("mu1_pt", &mu1_pt, &b_mu1_pt);
   fChain->SetBranchAddress("mu1_conept", &mu1_conept, &b_mu1_conept);
   fChain->SetBranchAddress("mu1_eta", &mu1_eta, &b_mu1_eta);
   fChain->SetBranchAddress("mu1_phi", &mu1_phi, &b_mu1_phi);
   fChain->SetBranchAddress("mu1_E", &mu1_E, &b_mu1_E);
   fChain->SetBranchAddress("mu1_charge", &mu1_charge, &b_mu1_charge);
   fChain->SetBranchAddress("mu1_miniRelIso", &mu1_miniRelIso, &b_mu1_miniRelIso);
   fChain->SetBranchAddress("mu1_PFRelIso04", &mu1_PFRelIso04, &b_mu1_PFRelIso04);
   fChain->SetBranchAddress("mu1_jetNDauChargedMVASel", &mu1_jetNDauChargedMVASel, &b_mu1_jetNDauChargedMVASel);
   fChain->SetBranchAddress("mu1_jetPtRel", &mu1_jetPtRel, &b_mu1_jetPtRel);
   fChain->SetBranchAddress("mu1_jetRelIso", &mu1_jetRelIso, &b_mu1_jetRelIso);
   fChain->SetBranchAddress("mu1_jetDeepJet", &mu1_jetDeepJet, &b_mu1_jetDeepJet);
   fChain->SetBranchAddress("mu1_sip3D", &mu1_sip3D, &b_mu1_sip3D);
   fChain->SetBranchAddress("mu1_dxy", &mu1_dxy, &b_mu1_dxy);
   fChain->SetBranchAddress("mu1_dxyAbs", &mu1_dxyAbs, &b_mu1_dxyAbs);
   fChain->SetBranchAddress("mu1_dz", &mu1_dz, &b_mu1_dz);
   fChain->SetBranchAddress("mu1_segmentCompatibility", &mu1_segmentCompatibility, &b_mu1_segmentCompatibility);
   fChain->SetBranchAddress("mu1_leptonMVA", &mu1_leptonMVA, &b_mu1_leptonMVA);
   fChain->SetBranchAddress("mu1_mediumID", &mu1_mediumID, &b_mu1_mediumID);
   fChain->SetBranchAddress("mu1_dpt_div_pt", &mu1_dpt_div_pt, &b_mu1_dpt_div_pt);
   fChain->SetBranchAddress("mu1_isfakeablesel", &mu1_isfakeablesel, &b_mu1_isfakeablesel);
   fChain->SetBranchAddress("mu1_ismvasel", &mu1_ismvasel, &b_mu1_ismvasel);
   fChain->SetBranchAddress("mu1_isGenMatched", &mu1_isGenMatched, &b_mu1_isGenMatched);
   fChain->SetBranchAddress("mu2_pt", &mu2_pt, &b_mu2_pt);
   fChain->SetBranchAddress("mu2_conept", &mu2_conept, &b_mu2_conept);
   fChain->SetBranchAddress("mu2_eta", &mu2_eta, &b_mu2_eta);
   fChain->SetBranchAddress("mu2_phi", &mu2_phi, &b_mu2_phi);
   fChain->SetBranchAddress("mu2_E", &mu2_E, &b_mu2_E);
   fChain->SetBranchAddress("mu2_charge", &mu2_charge, &b_mu2_charge);
   fChain->SetBranchAddress("mu2_miniRelIso", &mu2_miniRelIso, &b_mu2_miniRelIso);
   fChain->SetBranchAddress("mu2_PFRelIso04", &mu2_PFRelIso04, &b_mu2_PFRelIso04);
   fChain->SetBranchAddress("mu2_jetNDauChargedMVASel", &mu2_jetNDauChargedMVASel, &b_mu2_jetNDauChargedMVASel);
   fChain->SetBranchAddress("mu2_jetPtRel", &mu2_jetPtRel, &b_mu2_jetPtRel);
   fChain->SetBranchAddress("mu2_jetRelIso", &mu2_jetRelIso, &b_mu2_jetRelIso);
   fChain->SetBranchAddress("mu2_jetDeepJet", &mu2_jetDeepJet, &b_mu2_jetDeepJet);
   fChain->SetBranchAddress("mu2_sip3D", &mu2_sip3D, &b_mu2_sip3D);
   fChain->SetBranchAddress("mu2_dxy", &mu2_dxy, &b_mu2_dxy);
   fChain->SetBranchAddress("mu2_dxyAbs", &mu2_dxyAbs, &b_mu2_dxyAbs);
   fChain->SetBranchAddress("mu2_dz", &mu2_dz, &b_mu2_dz);
   fChain->SetBranchAddress("mu2_segmentCompatibility", &mu2_segmentCompatibility, &b_mu2_segmentCompatibility);
   fChain->SetBranchAddress("mu2_leptonMVA", &mu2_leptonMVA, &b_mu2_leptonMVA);
   fChain->SetBranchAddress("mu2_mediumID", &mu2_mediumID, &b_mu2_mediumID);
   fChain->SetBranchAddress("mu2_dpt_div_pt", &mu2_dpt_div_pt, &b_mu2_dpt_div_pt);
   fChain->SetBranchAddress("mu2_isfakeablesel", &mu2_isfakeablesel, &b_mu2_isfakeablesel);
   fChain->SetBranchAddress("mu2_ismvasel", &mu2_ismvasel, &b_mu2_ismvasel);
   fChain->SetBranchAddress("mu2_isGenMatched", &mu2_isGenMatched, &b_mu2_isGenMatched);
   fChain->SetBranchAddress("ele1_pt", &ele1_pt, &b_ele1_pt);
   fChain->SetBranchAddress("ele1_conept", &ele1_conept, &b_ele1_conept);
   fChain->SetBranchAddress("ele1_eta", &ele1_eta, &b_ele1_eta);
   fChain->SetBranchAddress("ele1_phi", &ele1_phi, &b_ele1_phi);
   fChain->SetBranchAddress("ele1_E", &ele1_E, &b_ele1_E);
   fChain->SetBranchAddress("ele1_charge", &ele1_charge, &b_ele1_charge);
   fChain->SetBranchAddress("ele1_miniRelIso", &ele1_miniRelIso, &b_ele1_miniRelIso);
   fChain->SetBranchAddress("ele1_PFRelIso04", &ele1_PFRelIso04, &b_ele1_PFRelIso04);
   fChain->SetBranchAddress("ele1_jetNDauChargedMVASel", &ele1_jetNDauChargedMVASel, &b_ele1_jetNDauChargedMVASel);
   fChain->SetBranchAddress("ele1_jetPtRel", &ele1_jetPtRel, &b_ele1_jetPtRel);
   fChain->SetBranchAddress("ele1_jetRelIso", &ele1_jetRelIso, &b_ele1_jetRelIso);
   fChain->SetBranchAddress("ele1_jetDeepJet", &ele1_jetDeepJet, &b_ele1_jetDeepJet);
   fChain->SetBranchAddress("ele1_sip3D", &ele1_sip3D, &b_ele1_sip3D);
   fChain->SetBranchAddress("ele1_dxy", &ele1_dxy, &b_ele1_dxy);
   fChain->SetBranchAddress("ele1_dxyAbs", &ele1_dxyAbs, &b_ele1_dxyAbs);
   fChain->SetBranchAddress("ele1_dz", &ele1_dz, &b_ele1_dz);
   fChain->SetBranchAddress("ele1_ntMVAeleID", &ele1_ntMVAeleID, &b_ele1_ntMVAeleID);
   fChain->SetBranchAddress("ele1_leptonMVA", &ele1_leptonMVA, &b_ele1_leptonMVA);
   fChain->SetBranchAddress("ele1_passesConversionVeto", &ele1_passesConversionVeto, &b_ele1_passesConversionVeto);
   fChain->SetBranchAddress("ele1_nMissingHits", &ele1_nMissingHits, &b_ele1_nMissingHits);
   fChain->SetBranchAddress("ele1_sigmaEtaEta", &ele1_sigmaEtaEta, &b_ele1_sigmaEtaEta);
   fChain->SetBranchAddress("ele1_HoE", &ele1_HoE, &b_ele1_HoE);
   fChain->SetBranchAddress("ele1_OoEminusOoP", &ele1_OoEminusOoP, &b_ele1_OoEminusOoP);
   fChain->SetBranchAddress("ele1_isfakeablesel", &ele1_isfakeablesel, &b_ele1_isfakeablesel);
   fChain->SetBranchAddress("ele1_ismvasel", &ele1_ismvasel, &b_ele1_ismvasel);
   fChain->SetBranchAddress("ele1_isGenMatched", &ele1_isGenMatched, &b_ele1_isGenMatched);
   fChain->SetBranchAddress("ele2_pt", &ele2_pt, &b_ele2_pt);
   fChain->SetBranchAddress("ele2_conept", &ele2_conept, &b_ele2_conept);
   fChain->SetBranchAddress("ele2_eta", &ele2_eta, &b_ele2_eta);
   fChain->SetBranchAddress("ele2_phi", &ele2_phi, &b_ele2_phi);
   fChain->SetBranchAddress("ele2_E", &ele2_E, &b_ele2_E);
   fChain->SetBranchAddress("ele2_charge", &ele2_charge, &b_ele2_charge);
   fChain->SetBranchAddress("ele2_miniRelIso", &ele2_miniRelIso, &b_ele2_miniRelIso);
   fChain->SetBranchAddress("ele2_PFRelIso04", &ele2_PFRelIso04, &b_ele2_PFRelIso04);
   fChain->SetBranchAddress("ele2_jetNDauChargedMVASel", &ele2_jetNDauChargedMVASel, &b_ele2_jetNDauChargedMVASel);
   fChain->SetBranchAddress("ele2_jetPtRel", &ele2_jetPtRel, &b_ele2_jetPtRel);
   fChain->SetBranchAddress("ele2_jetRelIso", &ele2_jetRelIso, &b_ele2_jetRelIso);
   fChain->SetBranchAddress("ele2_jetDeepJet", &ele2_jetDeepJet, &b_ele2_jetDeepJet);
   fChain->SetBranchAddress("ele2_sip3D", &ele2_sip3D, &b_ele2_sip3D);
   fChain->SetBranchAddress("ele2_dxy", &ele2_dxy, &b_ele2_dxy);
   fChain->SetBranchAddress("ele2_dxyAbs", &ele2_dxyAbs, &b_ele2_dxyAbs);
   fChain->SetBranchAddress("ele2_dz", &ele2_dz, &b_ele2_dz);
   fChain->SetBranchAddress("ele2_ntMVAeleID", &ele2_ntMVAeleID, &b_ele2_ntMVAeleID);
   fChain->SetBranchAddress("ele2_leptonMVA", &ele2_leptonMVA, &b_ele2_leptonMVA);
   fChain->SetBranchAddress("ele2_passesConversionVeto", &ele2_passesConversionVeto, &b_ele2_passesConversionVeto);
   fChain->SetBranchAddress("ele2_nMissingHits", &ele2_nMissingHits, &b_ele2_nMissingHits);
   fChain->SetBranchAddress("ele2_sigmaEtaEta", &ele2_sigmaEtaEta, &b_ele2_sigmaEtaEta);
   fChain->SetBranchAddress("ele2_HoE", &ele2_HoE, &b_ele2_HoE);
   fChain->SetBranchAddress("ele2_OoEminusOoP", &ele2_OoEminusOoP, &b_ele2_OoEminusOoP);
   fChain->SetBranchAddress("ele2_isfakeablesel", &ele2_isfakeablesel, &b_ele2_isfakeablesel);
   fChain->SetBranchAddress("ele2_ismvasel", &ele2_ismvasel, &b_ele2_ismvasel);
   fChain->SetBranchAddress("ele2_isGenMatched", &ele2_isGenMatched, &b_ele2_isGenMatched);
   fChain->SetBranchAddress("ak4Jet1_pt", &ak4Jet1_pt, &b_ak4Jet1_pt);
   fChain->SetBranchAddress("ak4Jet1_eta", &ak4Jet1_eta, &b_ak4Jet1_eta);
   fChain->SetBranchAddress("ak4Jet1_phi", &ak4Jet1_phi, &b_ak4Jet1_phi);
   fChain->SetBranchAddress("ak4Jet1_E", &ak4Jet1_E, &b_ak4Jet1_E);
   fChain->SetBranchAddress("ak4Jet1_CSV", &ak4Jet1_CSV, &b_ak4Jet1_CSV);
   fChain->SetBranchAddress("ak4Jet1_btagSF", &ak4Jet1_btagSF, &b_ak4Jet1_btagSF);
   fChain->SetBranchAddress("ak4Jet2_pt", &ak4Jet2_pt, &b_ak4Jet2_pt);
   fChain->SetBranchAddress("ak4Jet2_eta", &ak4Jet2_eta, &b_ak4Jet2_eta);
   fChain->SetBranchAddress("ak4Jet2_phi", &ak4Jet2_phi, &b_ak4Jet2_phi);
   fChain->SetBranchAddress("ak4Jet2_E", &ak4Jet2_E, &b_ak4Jet2_E);
   fChain->SetBranchAddress("ak4Jet2_CSV", &ak4Jet2_CSV, &b_ak4Jet2_CSV);
   fChain->SetBranchAddress("ak4Jet2_btagSF", &ak4Jet2_btagSF, &b_ak4Jet2_btagSF);
   fChain->SetBranchAddress("ak4Jet3_pt", &ak4Jet3_pt, &b_ak4Jet3_pt);
   fChain->SetBranchAddress("ak4Jet3_eta", &ak4Jet3_eta, &b_ak4Jet3_eta);
   fChain->SetBranchAddress("ak4Jet3_phi", &ak4Jet3_phi, &b_ak4Jet3_phi);
   fChain->SetBranchAddress("ak4Jet3_E", &ak4Jet3_E, &b_ak4Jet3_E);
   fChain->SetBranchAddress("ak4Jet3_CSV", &ak4Jet3_CSV, &b_ak4Jet3_CSV);
   fChain->SetBranchAddress("ak4Jet3_btagSF", &ak4Jet3_btagSF, &b_ak4Jet3_btagSF);
   fChain->SetBranchAddress("ak4Jet4_pt", &ak4Jet4_pt, &b_ak4Jet4_pt);
   fChain->SetBranchAddress("ak4Jet4_eta", &ak4Jet4_eta, &b_ak4Jet4_eta);
   fChain->SetBranchAddress("ak4Jet4_phi", &ak4Jet4_phi, &b_ak4Jet4_phi);
   fChain->SetBranchAddress("ak4Jet4_E", &ak4Jet4_E, &b_ak4Jet4_E);
   fChain->SetBranchAddress("ak4Jet4_CSV", &ak4Jet4_CSV, &b_ak4Jet4_CSV);
   fChain->SetBranchAddress("ak4Jet4_btagSF", &ak4Jet4_btagSF, &b_ak4Jet4_btagSF);
   fChain->SetBranchAddress("ak4JetVBF1_pt", &ak4JetVBF1_pt, &b_ak4JetVBF1_pt);
   fChain->SetBranchAddress("ak4JetVBF1_eta", &ak4JetVBF1_eta, &b_ak4JetVBF1_eta);
   fChain->SetBranchAddress("ak4JetVBF1_phi", &ak4JetVBF1_phi, &b_ak4JetVBF1_phi);
   fChain->SetBranchAddress("ak4JetVBF1_E", &ak4JetVBF1_E, &b_ak4JetVBF1_E);
   fChain->SetBranchAddress("ak4JetVBF1_CSV", &ak4JetVBF1_CSV, &b_ak4JetVBF1_CSV);
   fChain->SetBranchAddress("ak4JetVBF1_btagSF", &ak4JetVBF1_btagSF, &b_ak4JetVBF1_btagSF);
   fChain->SetBranchAddress("ak4JetVBF2_pt", &ak4JetVBF2_pt, &b_ak4JetVBF2_pt);
   fChain->SetBranchAddress("ak4JetVBF2_eta", &ak4JetVBF2_eta, &b_ak4JetVBF2_eta);
   fChain->SetBranchAddress("ak4JetVBF2_phi", &ak4JetVBF2_phi, &b_ak4JetVBF2_phi);
   fChain->SetBranchAddress("ak4JetVBF2_E", &ak4JetVBF2_E, &b_ak4JetVBF2_E);
   fChain->SetBranchAddress("ak4JetVBF2_CSV", &ak4JetVBF2_CSV, &b_ak4JetVBF2_CSV);
   fChain->SetBranchAddress("ak4JetVBF2_btagSF", &ak4JetVBF2_btagSF, &b_ak4JetVBF2_btagSF);
   fChain->SetBranchAddress("ak8Jet1_pt", &ak8Jet1_pt, &b_ak8Jet1_pt);
   fChain->SetBranchAddress("ak8Jet1_eta", &ak8Jet1_eta, &b_ak8Jet1_eta);
   fChain->SetBranchAddress("ak8Jet1_phi", &ak8Jet1_phi, &b_ak8Jet1_phi);
   fChain->SetBranchAddress("ak8Jet1_E", &ak8Jet1_E, &b_ak8Jet1_E);
   fChain->SetBranchAddress("ak8Jet1_msoftdrop", &ak8Jet1_msoftdrop, &b_ak8Jet1_msoftdrop);
   fChain->SetBranchAddress("ak8Jet1_tau1", &ak8Jet1_tau1, &b_ak8Jet1_tau1);
   fChain->SetBranchAddress("ak8Jet1_tau2", &ak8Jet1_tau2, &b_ak8Jet1_tau2);
   fChain->SetBranchAddress("ak8Jet1_subjet0_pt", &ak8Jet1_subjet0_pt, &b_ak8Jet1_subjet0_pt);
   fChain->SetBranchAddress("ak8Jet1_subjet0_eta", &ak8Jet1_subjet0_eta, &b_ak8Jet1_subjet0_eta);
   fChain->SetBranchAddress("ak8Jet1_subjet0_phi", &ak8Jet1_subjet0_phi, &b_ak8Jet1_subjet0_phi);
   fChain->SetBranchAddress("ak8Jet1_subjet0_CSV", &ak8Jet1_subjet0_CSV, &b_ak8Jet1_subjet0_CSV);
   fChain->SetBranchAddress("ak8Jet1_subjet1_pt", &ak8Jet1_subjet1_pt, &b_ak8Jet1_subjet1_pt);
   fChain->SetBranchAddress("ak8Jet1_subjet1_eta", &ak8Jet1_subjet1_eta, &b_ak8Jet1_subjet1_eta);
   fChain->SetBranchAddress("ak8Jet1_subjet1_phi", &ak8Jet1_subjet1_phi, &b_ak8Jet1_subjet1_phi);
   fChain->SetBranchAddress("ak8Jet1_subjet1_CSV", &ak8Jet1_subjet1_CSV, &b_ak8Jet1_subjet1_CSV);
   fChain->SetBranchAddress("ak8Jet2_pt", &ak8Jet2_pt, &b_ak8Jet2_pt);
   fChain->SetBranchAddress("ak8Jet2_eta", &ak8Jet2_eta, &b_ak8Jet2_eta);
   fChain->SetBranchAddress("ak8Jet2_phi", &ak8Jet2_phi, &b_ak8Jet2_phi);
   fChain->SetBranchAddress("ak8Jet2_E", &ak8Jet2_E, &b_ak8Jet2_E);
   fChain->SetBranchAddress("ak8Jet2_msoftdrop", &ak8Jet2_msoftdrop, &b_ak8Jet2_msoftdrop);
   fChain->SetBranchAddress("ak8Jet2_tau1", &ak8Jet2_tau1, &b_ak8Jet2_tau1);
   fChain->SetBranchAddress("ak8Jet2_tau2", &ak8Jet2_tau2, &b_ak8Jet2_tau2);
   fChain->SetBranchAddress("ak8Jet2_subjet0_pt", &ak8Jet2_subjet0_pt, &b_ak8Jet2_subjet0_pt);
   fChain->SetBranchAddress("ak8Jet2_subjet0_eta", &ak8Jet2_subjet0_eta, &b_ak8Jet2_subjet0_eta);
   fChain->SetBranchAddress("ak8Jet2_subjet0_phi", &ak8Jet2_subjet0_phi, &b_ak8Jet2_subjet0_phi);
   fChain->SetBranchAddress("ak8Jet2_subjet0_CSV", &ak8Jet2_subjet0_CSV, &b_ak8Jet2_subjet0_CSV);
   fChain->SetBranchAddress("ak8Jet2_subjet1_pt", &ak8Jet2_subjet1_pt, &b_ak8Jet2_subjet1_pt);
   fChain->SetBranchAddress("ak8Jet2_subjet1_eta", &ak8Jet2_subjet1_eta, &b_ak8Jet2_subjet1_eta);
   fChain->SetBranchAddress("ak8Jet2_subjet1_phi", &ak8Jet2_subjet1_phi, &b_ak8Jet2_subjet1_phi);
   fChain->SetBranchAddress("ak8Jet2_subjet1_CSV", &ak8Jet2_subjet1_CSV, &b_ak8Jet2_subjet1_CSV);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETphi", &PFMETphi, &b_PFMETphi);
   fChain->SetBranchAddress("HME", &HME, &b_HME);
   fChain->SetBranchAddress("PU_weight", &PU_weight, &b_PU_weight);
   fChain->SetBranchAddress("PU_jetID_SF", &PU_jetID_SF, &b_PU_jetID_SF);
   fChain->SetBranchAddress("MC_weight", &MC_weight, &b_MC_weight);
   fChain->SetBranchAddress("topPt_wgt", &topPt_wgt, &b_topPt_wgt);
   fChain->SetBranchAddress("btag_SF", &btag_SF, &b_btag_SF);
   fChain->SetBranchAddress("trigger_SF", &trigger_SF, &b_trigger_SF);
   fChain->SetBranchAddress("lepton_IDSF", &lepton_IDSF, &b_lepton_IDSF);
   fChain->SetBranchAddress("lepton_IDSF_recoToLoose", &lepton_IDSF_recoToLoose, &b_lepton_IDSF_recoToLoose);
   fChain->SetBranchAddress("lepton_IDSF_looseToTight", &lepton_IDSF_looseToTight, &b_lepton_IDSF_looseToTight);
   fChain->SetBranchAddress("L1prefire", &L1prefire, &b_L1prefire);
   fChain->SetBranchAddress("fakeRate", &fakeRate, &b_fakeRate);
   fChain->SetBranchAddress("vbf_m_jj", &vbf_m_jj, &b_vbf_m_jj);
   fChain->SetBranchAddress("vbf_dEta_jj", &vbf_dEta_jj, &b_vbf_dEta_jj);
   fChain->SetBranchAddress("trigger_SF_down", &trigger_SF_down, &b_trigger_SF_down);
   fChain->SetBranchAddress("trigger_SF_up", &trigger_SF_up, &b_trigger_SF_up);
   fChain->SetBranchAddress("lepton1_{branchname}_nominal", &lepton1_{branchname}_nominal, &b_lepton1_{branchname}_nominal);
   fChain->SetBranchAddress("lepton1_{branchname}_up", &lepton1_{branchname}_up, &b_lepton1_{branchname}_up);
   fChain->SetBranchAddress("lepton1_{branchname}_down", &lepton1_{branchname}_down, &b_lepton1_{branchname}_down);
   fChain->SetBranchAddress("lepton2_{branchname}_nominal", &lepton2_{branchname}_nominal, &b_lepton2_{branchname}_nominal);
   fChain->SetBranchAddress("lepton2_{branchname}_up", &lepton2_{branchname}_up, &b_lepton2_{branchname}_up);
   fChain->SetBranchAddress("lepton2_{branchname}_down", &lepton2_{branchname}_down, &b_lepton2_{branchname}_down);
   fChain->SetBranchAddress("genMET_pT", &genMET_pT, &b_genMET_pT);
   fChain->SetBranchAddress("genMET_phi", &genMET_phi, &b_genMET_phi);
   fChain->SetBranchAddress("met_covxx", &met_covxx, &b_met_covxx);
   fChain->SetBranchAddress("met_covyy", &met_covyy, &b_met_covyy);
   fChain->SetBranchAddress("met_covxy", &met_covxy, &b_met_covxy);
   fChain->SetBranchAddress("findAllGen", &findAllGen, &b_findAllGen);
   fChain->SetBranchAddress("finalxtohhtobbww", &finalxtohhtobbww, &b_finalxtohhtobbww);
   fChain->SetBranchAddress("is_genboosted", &is_genboosted, &b_is_genboosted);
   fChain->SetBranchAddress("is_genDL", &is_genDL, &b_is_genDL);
   fChain->SetBranchAddress("is_genSL", &is_genSL, &b_is_genSL);
   fChain->SetBranchAddress("is_genee", &is_genee, &b_is_genee);
   fChain->SetBranchAddress("is_genem", &is_genem, &b_is_genem);
   fChain->SetBranchAddress("is_genmm", &is_genmm, &b_is_genmm);
   fChain->SetBranchAddress("genhtobb_pt", &genhtobb_pt, &b_genhtobb_pt);
   fChain->SetBranchAddress("genhtobb_eta", &genhtobb_eta, &b_genhtobb_eta);
   fChain->SetBranchAddress("genhtobb_phi", &genhtobb_phi, &b_genhtobb_phi);
   fChain->SetBranchAddress("genhtobb_mass", &genhtobb_mass, &b_genhtobb_mass);
   fChain->SetBranchAddress("genhtobb_index", &genhtobb_index, &b_genhtobb_index);
   fChain->SetBranchAddress("genhtobb_gendR", &genhtobb_gendR, &b_genhtobb_gendR);
   fChain->SetBranchAddress("genhtobb_matchidx", &genhtobb_matchidx, &b_genhtobb_matchidx);
   fChain->SetBranchAddress("genhtobb_genPartIdxMother", &genhtobb_genPartIdxMother, &b_genhtobb_genPartIdxMother);
   fChain->SetBranchAddress("genhtobb_pdgId", &genhtobb_pdgId, &b_genhtobb_pdgId);
   fChain->SetBranchAddress("genhtobb_statusFlags", &genhtobb_statusFlags, &b_genhtobb_statusFlags);
   fChain->SetBranchAddress("genbjet1_pt", &genbjet1_pt, &b_genbjet1_pt);
   fChain->SetBranchAddress("genbjet1_eta", &genbjet1_eta, &b_genbjet1_eta);
   fChain->SetBranchAddress("genbjet1_phi", &genbjet1_phi, &b_genbjet1_phi);
   fChain->SetBranchAddress("genbjet1_mass", &genbjet1_mass, &b_genbjet1_mass);
   fChain->SetBranchAddress("genbjet1_index", &genbjet1_index, &b_genbjet1_index);
   fChain->SetBranchAddress("genbjet1_gendR", &genbjet1_gendR, &b_genbjet1_gendR);
   fChain->SetBranchAddress("genbjet1_matchidx", &genbjet1_matchidx, &b_genbjet1_matchidx);
   fChain->SetBranchAddress("genbjet1_hadronFlavour", &genbjet1_hadronFlavour, &b_genbjet1_hadronFlavour);
   fChain->SetBranchAddress("genbjet1_partonFlavour", &genbjet1_partonFlavour, &b_genbjet1_partonFlavour);
   fChain->SetBranchAddress("genbjet1_partondR", &genbjet1_partondR, &b_genbjet1_partondR);
   fChain->SetBranchAddress("genq2_pt", &genq2_pt, &b_genq2_pt);
   fChain->SetBranchAddress("genq2_eta", &genq2_eta, &b_genq2_eta);
   fChain->SetBranchAddress("genq2_phi", &genq2_phi, &b_genq2_phi);
   fChain->SetBranchAddress("genq2_mass", &genq2_mass, &b_genq2_mass);
   fChain->SetBranchAddress("genq2_index", &genq2_index, &b_genq2_index);
   fChain->SetBranchAddress("genq2_gendR", &genq2_gendR, &b_genq2_gendR);
   fChain->SetBranchAddress("genq2_matchidx", &genq2_matchidx, &b_genq2_matchidx);
   fChain->SetBranchAddress("genq2_genPartIdxMother", &genq2_genPartIdxMother, &b_genq2_genPartIdxMother);
   fChain->SetBranchAddress("genq2_pdgId", &genq2_pdgId, &b_genq2_pdgId);
   fChain->SetBranchAddress("genq2_statusFlags", &genq2_statusFlags, &b_genq2_statusFlags);
   fChain->SetBranchAddress("genbjet2nu_pt", &genbjet2nu_pt, &b_genbjet2nu_pt);
   fChain->SetBranchAddress("genbjet2nu_eta", &genbjet2nu_eta, &b_genbjet2nu_eta);
   fChain->SetBranchAddress("genbjet2nu_phi", &genbjet2nu_phi, &b_genbjet2nu_phi);
   fChain->SetBranchAddress("genbjet2nu_mass", &genbjet2nu_mass, &b_genbjet2nu_mass);
   fChain->SetBranchAddress("genbjet2nu_index", &genbjet2nu_index, &b_genbjet2nu_index);
   fChain->SetBranchAddress("genbjet2nu_gendR", &genbjet2nu_gendR, &b_genbjet2nu_gendR);
   fChain->SetBranchAddress("genbjet2nu_matchidx", &genbjet2nu_matchidx, &b_genbjet2nu_matchidx);
   fChain->SetBranchAddress("genbjet2nu_genPartIdxMother", &genbjet2nu_genPartIdxMother, &b_genbjet2nu_genPartIdxMother);
   fChain->SetBranchAddress("genbjet2nu_pdgId", &genbjet2nu_pdgId, &b_genbjet2nu_pdgId);
   fChain->SetBranchAddress("genbjet2nu_statusFlags", &genbjet2nu_statusFlags, &b_genbjet2nu_statusFlags);
   fChain->SetBranchAddress("genak8jet_pt", &genak8jet_pt, &b_genak8jet_pt);
   fChain->SetBranchAddress("genak8jet_eta", &genak8jet_eta, &b_genak8jet_eta);
   fChain->SetBranchAddress("genak8jet_phi", &genak8jet_phi, &b_genak8jet_phi);
   fChain->SetBranchAddress("genak8jet_mass", &genak8jet_mass, &b_genak8jet_mass);
   fChain->SetBranchAddress("genak8jet_index", &genak8jet_index, &b_genak8jet_index);
   fChain->SetBranchAddress("genak8jet_gendR", &genak8jet_gendR, &b_genak8jet_gendR);
   fChain->SetBranchAddress("genak8jet_matchidx", &genak8jet_matchidx, &b_genak8jet_matchidx);
   fChain->SetBranchAddress("genak8jet_hadronFlavour", &genak8jet_hadronFlavour, &b_genak8jet_hadronFlavour);
   fChain->SetBranchAddress("genak8jet_partonFlavour", &genak8jet_partonFlavour, &b_genak8jet_partonFlavour);
   fChain->SetBranchAddress("genak8jet_partondR", &genak8jet_partondR, &b_genak8jet_partondR);
   fChain->SetBranchAddress("genbjet1nu_pt", &genbjet1nu_pt, &b_genbjet1nu_pt);
   fChain->SetBranchAddress("genbjet1nu_eta", &genbjet1nu_eta, &b_genbjet1nu_eta);
   fChain->SetBranchAddress("genbjet1nu_phi", &genbjet1nu_phi, &b_genbjet1nu_phi);
   fChain->SetBranchAddress("genbjet1nu_mass", &genbjet1nu_mass, &b_genbjet1nu_mass);
   fChain->SetBranchAddress("genbjet1nu_index", &genbjet1nu_index, &b_genbjet1nu_index);
   fChain->SetBranchAddress("genbjet1nu_gendR", &genbjet1nu_gendR, &b_genbjet1nu_gendR);
   fChain->SetBranchAddress("genbjet1nu_matchidx", &genbjet1nu_matchidx, &b_genbjet1nu_matchidx);
   fChain->SetBranchAddress("genbjet1nu_genPartIdxMother", &genbjet1nu_genPartIdxMother, &b_genbjet1nu_genPartIdxMother);
   fChain->SetBranchAddress("genbjet1nu_pdgId", &genbjet1nu_pdgId, &b_genbjet1nu_pdgId);
   fChain->SetBranchAddress("genbjet1nu_statusFlags", &genbjet1nu_statusFlags, &b_genbjet1nu_statusFlags);
   fChain->SetBranchAddress("gennu1_pt", &gennu1_pt, &b_gennu1_pt);
   fChain->SetBranchAddress("gennu1_eta", &gennu1_eta, &b_gennu1_eta);
   fChain->SetBranchAddress("gennu1_phi", &gennu1_phi, &b_gennu1_phi);
   fChain->SetBranchAddress("gennu1_mass", &gennu1_mass, &b_gennu1_mass);
   fChain->SetBranchAddress("gennu1_index", &gennu1_index, &b_gennu1_index);
   fChain->SetBranchAddress("gennu1_gendR", &gennu1_gendR, &b_gennu1_gendR);
   fChain->SetBranchAddress("gennu1_matchidx", &gennu1_matchidx, &b_gennu1_matchidx);
   fChain->SetBranchAddress("gennu1_genPartIdxMother", &gennu1_genPartIdxMother, &b_gennu1_genPartIdxMother);
   fChain->SetBranchAddress("gennu1_pdgId", &gennu1_pdgId, &b_gennu1_pdgId);
   fChain->SetBranchAddress("gennu1_statusFlags", &gennu1_statusFlags, &b_gennu1_statusFlags);
   fChain->SetBranchAddress("genbjet2_pt", &genbjet2_pt, &b_genbjet2_pt);
   fChain->SetBranchAddress("genbjet2_eta", &genbjet2_eta, &b_genbjet2_eta);
   fChain->SetBranchAddress("genbjet2_phi", &genbjet2_phi, &b_genbjet2_phi);
   fChain->SetBranchAddress("genbjet2_mass", &genbjet2_mass, &b_genbjet2_mass);
   fChain->SetBranchAddress("genbjet2_index", &genbjet2_index, &b_genbjet2_index);
   fChain->SetBranchAddress("genbjet2_gendR", &genbjet2_gendR, &b_genbjet2_gendR);
   fChain->SetBranchAddress("genbjet2_matchidx", &genbjet2_matchidx, &b_genbjet2_matchidx);
   fChain->SetBranchAddress("genbjet2_hadronFlavour", &genbjet2_hadronFlavour, &b_genbjet2_hadronFlavour);
   fChain->SetBranchAddress("genbjet2_partonFlavour", &genbjet2_partonFlavour, &b_genbjet2_partonFlavour);
   fChain->SetBranchAddress("genbjet2_partondR", &genbjet2_partondR, &b_genbjet2_partondR);
   fChain->SetBranchAddress("genb2_pt", &genb2_pt, &b_genb2_pt);
   fChain->SetBranchAddress("genb2_eta", &genb2_eta, &b_genb2_eta);
   fChain->SetBranchAddress("genb2_phi", &genb2_phi, &b_genb2_phi);
   fChain->SetBranchAddress("genb2_mass", &genb2_mass, &b_genb2_mass);
   fChain->SetBranchAddress("genb2_index", &genb2_index, &b_genb2_index);
   fChain->SetBranchAddress("genb2_gendR", &genb2_gendR, &b_genb2_gendR);
   fChain->SetBranchAddress("genb2_matchidx", &genb2_matchidx, &b_genb2_matchidx);
   fChain->SetBranchAddress("genb2_genPartIdxMother", &genb2_genPartIdxMother, &b_genb2_genPartIdxMother);
   fChain->SetBranchAddress("genb2_pdgId", &genb2_pdgId, &b_genb2_pdgId);
   fChain->SetBranchAddress("genb2_statusFlags", &genb2_statusFlags, &b_genb2_statusFlags);
   fChain->SetBranchAddress("gennu2_pt", &gennu2_pt, &b_gennu2_pt);
   fChain->SetBranchAddress("gennu2_eta", &gennu2_eta, &b_gennu2_eta);
   fChain->SetBranchAddress("gennu2_phi", &gennu2_phi, &b_gennu2_phi);
   fChain->SetBranchAddress("gennu2_mass", &gennu2_mass, &b_gennu2_mass);
   fChain->SetBranchAddress("gennu2_index", &gennu2_index, &b_gennu2_index);
   fChain->SetBranchAddress("gennu2_gendR", &gennu2_gendR, &b_gennu2_gendR);
   fChain->SetBranchAddress("gennu2_matchidx", &gennu2_matchidx, &b_gennu2_matchidx);
   fChain->SetBranchAddress("gennu2_genPartIdxMother", &gennu2_genPartIdxMother, &b_gennu2_genPartIdxMother);
   fChain->SetBranchAddress("gennu2_pdgId", &gennu2_pdgId, &b_gennu2_pdgId);
   fChain->SetBranchAddress("gennu2_statusFlags", &gennu2_statusFlags, &b_gennu2_statusFlags);
   fChain->SetBranchAddress("genw2_pt", &genw2_pt, &b_genw2_pt);
   fChain->SetBranchAddress("genw2_eta", &genw2_eta, &b_genw2_eta);
   fChain->SetBranchAddress("genw2_phi", &genw2_phi, &b_genw2_phi);
   fChain->SetBranchAddress("genw2_mass", &genw2_mass, &b_genw2_mass);
   fChain->SetBranchAddress("genw2_index", &genw2_index, &b_genw2_index);
   fChain->SetBranchAddress("genw2_gendR", &genw2_gendR, &b_genw2_gendR);
   fChain->SetBranchAddress("genw2_matchidx", &genw2_matchidx, &b_genw2_matchidx);
   fChain->SetBranchAddress("genw2_genPartIdxMother", &genw2_genPartIdxMother, &b_genw2_genPartIdxMother);
   fChain->SetBranchAddress("genw2_pdgId", &genw2_pdgId, &b_genw2_pdgId);
   fChain->SetBranchAddress("genw2_statusFlags", &genw2_statusFlags, &b_genw2_statusFlags);
   fChain->SetBranchAddress("genqjet1_pt", &genqjet1_pt, &b_genqjet1_pt);
   fChain->SetBranchAddress("genqjet1_eta", &genqjet1_eta, &b_genqjet1_eta);
   fChain->SetBranchAddress("genqjet1_phi", &genqjet1_phi, &b_genqjet1_phi);
   fChain->SetBranchAddress("genqjet1_mass", &genqjet1_mass, &b_genqjet1_mass);
   fChain->SetBranchAddress("genqjet1_index", &genqjet1_index, &b_genqjet1_index);
   fChain->SetBranchAddress("genqjet1_gendR", &genqjet1_gendR, &b_genqjet1_gendR);
   fChain->SetBranchAddress("genqjet1_matchidx", &genqjet1_matchidx, &b_genqjet1_matchidx);
   fChain->SetBranchAddress("genqjet1_hadronFlavour", &genqjet1_hadronFlavour, &b_genqjet1_hadronFlavour);
   fChain->SetBranchAddress("genqjet1_partonFlavour", &genqjet1_partonFlavour, &b_genqjet1_partonFlavour);
   fChain->SetBranchAddress("genqjet1_partondR", &genqjet1_partondR, &b_genqjet1_partondR);
   fChain->SetBranchAddress("genqjet2_pt", &genqjet2_pt, &b_genqjet2_pt);
   fChain->SetBranchAddress("genqjet2_eta", &genqjet2_eta, &b_genqjet2_eta);
   fChain->SetBranchAddress("genqjet2_phi", &genqjet2_phi, &b_genqjet2_phi);
   fChain->SetBranchAddress("genqjet2_mass", &genqjet2_mass, &b_genqjet2_mass);
   fChain->SetBranchAddress("genqjet2_index", &genqjet2_index, &b_genqjet2_index);
   fChain->SetBranchAddress("genqjet2_gendR", &genqjet2_gendR, &b_genqjet2_gendR);
   fChain->SetBranchAddress("genqjet2_matchidx", &genqjet2_matchidx, &b_genqjet2_matchidx);
   fChain->SetBranchAddress("genqjet2_hadronFlavour", &genqjet2_hadronFlavour, &b_genqjet2_hadronFlavour);
   fChain->SetBranchAddress("genqjet2_partonFlavour", &genqjet2_partonFlavour, &b_genqjet2_partonFlavour);
   fChain->SetBranchAddress("genqjet2_partondR", &genqjet2_partondR, &b_genqjet2_partondR);
   fChain->SetBranchAddress("genw1_pt", &genw1_pt, &b_genw1_pt);
   fChain->SetBranchAddress("genw1_eta", &genw1_eta, &b_genw1_eta);
   fChain->SetBranchAddress("genw1_phi", &genw1_phi, &b_genw1_phi);
   fChain->SetBranchAddress("genw1_mass", &genw1_mass, &b_genw1_mass);
   fChain->SetBranchAddress("genw1_index", &genw1_index, &b_genw1_index);
   fChain->SetBranchAddress("genw1_gendR", &genw1_gendR, &b_genw1_gendR);
   fChain->SetBranchAddress("genw1_matchidx", &genw1_matchidx, &b_genw1_matchidx);
   fChain->SetBranchAddress("genw1_genPartIdxMother", &genw1_genPartIdxMother, &b_genw1_genPartIdxMother);
   fChain->SetBranchAddress("genw1_pdgId", &genw1_pdgId, &b_genw1_pdgId);
   fChain->SetBranchAddress("genw1_statusFlags", &genw1_statusFlags, &b_genw1_statusFlags);
   fChain->SetBranchAddress("genl2_pt", &genl2_pt, &b_genl2_pt);
   fChain->SetBranchAddress("genl2_eta", &genl2_eta, &b_genl2_eta);
   fChain->SetBranchAddress("genl2_phi", &genl2_phi, &b_genl2_phi);
   fChain->SetBranchAddress("genl2_mass", &genl2_mass, &b_genl2_mass);
   fChain->SetBranchAddress("genl2_index", &genl2_index, &b_genl2_index);
   fChain->SetBranchAddress("genl2_gendR", &genl2_gendR, &b_genl2_gendR);
   fChain->SetBranchAddress("genl2_matchidx", &genl2_matchidx, &b_genl2_matchidx);
   fChain->SetBranchAddress("genl2_genPartIdxMother", &genl2_genPartIdxMother, &b_genl2_genPartIdxMother);
   fChain->SetBranchAddress("genl2_pdgId", &genl2_pdgId, &b_genl2_pdgId);
   fChain->SetBranchAddress("genl2_statusFlags", &genl2_statusFlags, &b_genl2_statusFlags);
   fChain->SetBranchAddress("genl1_pt", &genl1_pt, &b_genl1_pt);
   fChain->SetBranchAddress("genl1_eta", &genl1_eta, &b_genl1_eta);
   fChain->SetBranchAddress("genl1_phi", &genl1_phi, &b_genl1_phi);
   fChain->SetBranchAddress("genl1_mass", &genl1_mass, &b_genl1_mass);
   fChain->SetBranchAddress("genl1_index", &genl1_index, &b_genl1_index);
   fChain->SetBranchAddress("genl1_gendR", &genl1_gendR, &b_genl1_gendR);
   fChain->SetBranchAddress("genl1_matchidx", &genl1_matchidx, &b_genl1_matchidx);
   fChain->SetBranchAddress("genl1_genPartIdxMother", &genl1_genPartIdxMother, &b_genl1_genPartIdxMother);
   fChain->SetBranchAddress("genl1_pdgId", &genl1_pdgId, &b_genl1_pdgId);
   fChain->SetBranchAddress("genl1_statusFlags", &genl1_statusFlags, &b_genl1_statusFlags);
   fChain->SetBranchAddress("genxtohh_pt", &genxtohh_pt, &b_genxtohh_pt);
   fChain->SetBranchAddress("genxtohh_eta", &genxtohh_eta, &b_genxtohh_eta);
   fChain->SetBranchAddress("genxtohh_phi", &genxtohh_phi, &b_genxtohh_phi);
   fChain->SetBranchAddress("genxtohh_mass", &genxtohh_mass, &b_genxtohh_mass);
   fChain->SetBranchAddress("genxtohh_index", &genxtohh_index, &b_genxtohh_index);
   fChain->SetBranchAddress("genxtohh_gendR", &genxtohh_gendR, &b_genxtohh_gendR);
   fChain->SetBranchAddress("genxtohh_matchidx", &genxtohh_matchidx, &b_genxtohh_matchidx);
   fChain->SetBranchAddress("genxtohh_genPartIdxMother", &genxtohh_genPartIdxMother, &b_genxtohh_genPartIdxMother);
   fChain->SetBranchAddress("genxtohh_pdgId", &genxtohh_pdgId, &b_genxtohh_pdgId);
   fChain->SetBranchAddress("genxtohh_statusFlags", &genxtohh_statusFlags, &b_genxtohh_statusFlags);
   fChain->SetBranchAddress("genhtoww_pt", &genhtoww_pt, &b_genhtoww_pt);
   fChain->SetBranchAddress("genhtoww_eta", &genhtoww_eta, &b_genhtoww_eta);
   fChain->SetBranchAddress("genhtoww_phi", &genhtoww_phi, &b_genhtoww_phi);
   fChain->SetBranchAddress("genhtoww_mass", &genhtoww_mass, &b_genhtoww_mass);
   fChain->SetBranchAddress("genhtoww_index", &genhtoww_index, &b_genhtoww_index);
   fChain->SetBranchAddress("genhtoww_gendR", &genhtoww_gendR, &b_genhtoww_gendR);
   fChain->SetBranchAddress("genhtoww_matchidx", &genhtoww_matchidx, &b_genhtoww_matchidx);
   fChain->SetBranchAddress("genhtoww_genPartIdxMother", &genhtoww_genPartIdxMother, &b_genhtoww_genPartIdxMother);
   fChain->SetBranchAddress("genhtoww_pdgId", &genhtoww_pdgId, &b_genhtoww_pdgId);
   fChain->SetBranchAddress("genhtoww_statusFlags", &genhtoww_statusFlags, &b_genhtoww_statusFlags);
   fChain->SetBranchAddress("genb1_pt", &genb1_pt, &b_genb1_pt);
   fChain->SetBranchAddress("genb1_eta", &genb1_eta, &b_genb1_eta);
   fChain->SetBranchAddress("genb1_phi", &genb1_phi, &b_genb1_phi);
   fChain->SetBranchAddress("genb1_mass", &genb1_mass, &b_genb1_mass);
   fChain->SetBranchAddress("genb1_index", &genb1_index, &b_genb1_index);
   fChain->SetBranchAddress("genb1_gendR", &genb1_gendR, &b_genb1_gendR);
   fChain->SetBranchAddress("genb1_matchidx", &genb1_matchidx, &b_genb1_matchidx);
   fChain->SetBranchAddress("genb1_genPartIdxMother", &genb1_genPartIdxMother, &b_genb1_genPartIdxMother);
   fChain->SetBranchAddress("genb1_pdgId", &genb1_pdgId, &b_genb1_pdgId);
   fChain->SetBranchAddress("genb1_statusFlags", &genb1_statusFlags, &b_genb1_statusFlags);
   fChain->SetBranchAddress("genq1_pt", &genq1_pt, &b_genq1_pt);
   fChain->SetBranchAddress("genq1_eta", &genq1_eta, &b_genq1_eta);
   fChain->SetBranchAddress("genq1_phi", &genq1_phi, &b_genq1_phi);
   fChain->SetBranchAddress("genq1_mass", &genq1_mass, &b_genq1_mass);
   fChain->SetBranchAddress("genq1_index", &genq1_index, &b_genq1_index);
   fChain->SetBranchAddress("genq1_gendR", &genq1_gendR, &b_genq1_gendR);
   fChain->SetBranchAddress("genq1_matchidx", &genq1_matchidx, &b_genq1_matchidx);
   fChain->SetBranchAddress("genq1_genPartIdxMother", &genq1_genPartIdxMother, &b_genq1_genPartIdxMother);
   fChain->SetBranchAddress("genq1_pdgId", &genq1_pdgId, &b_genq1_pdgId);
   fChain->SetBranchAddress("genq1_statusFlags", &genq1_statusFlags, &b_genq1_statusFlags);
   Notify();
}

Bool_t syncTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void syncTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t syncTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef syncTree_cxx
