#ifndef EVENT_HPP
#define EVENT_HPP

#include <array>

#include "TTree.h"
#include "ReadoutConstants.hpp"

struct Event
{   
    void ConnectTree(TTree* tree);

    //reco objects
    // 8x4x20 = 640
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_pt = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_eta = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_phi = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_mass = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_corr = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_res = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_btag = {};
    std::array<Float_t, ReadConst::MAX_RECO_JET> reco_jet_qvg = {};

    // 6x4x2 = 48
    std::array<Float_t, ReadConst::MAX_RECO_LEP> reco_lep_pt = {};
    std::array<Float_t, ReadConst::MAX_RECO_LEP> reco_lep_eta = {};
    std::array<Float_t, ReadConst::MAX_RECO_LEP> reco_lep_phi = {};
    std::array<Float_t, ReadConst::MAX_RECO_LEP> reco_lep_mass = {};
    std::array<Int_t, ReadConst::MAX_RECO_LEP> reco_lep_type = {};
    std::array<Int_t, ReadConst::MAX_RECO_LEP> reco_lep_gen_kind = {};

    // reco fatjets
    // 5x4x5 = 100
    std::array<Float_t, ReadConst::MAX_RECO_FAT_JET> reco_fatjet_pt = {};
    std::array<Float_t, ReadConst::MAX_RECO_FAT_JET> reco_fatjet_eta = {};
    std::array<Float_t, ReadConst::MAX_RECO_FAT_JET> reco_fatjet_phi = {};
    std::array<Float_t, ReadConst::MAX_RECO_FAT_JET> reco_fatjet_mass = {};
    std::array<Float_t, ReadConst::MAX_RECO_FAT_JET> reco_fatjet_btag_HbbvsQCD = {};

    // gen objects
    // 4x4x4 = 64
    std::array<Float_t, ReadConst::MAX_GEN_QUARK> gen_quark_pt = {};
    std::array<Float_t, ReadConst::MAX_GEN_QUARK> gen_quark_eta = {};
    std::array<Float_t, ReadConst::MAX_GEN_QUARK> gen_quark_phi = {};
    std::array<Float_t, ReadConst::MAX_GEN_QUARK> gen_quark_mass = {};

    // 4x4x4 = 64
    std::array<Float_t, ReadConst::MAX_GEN_LEP> gen_lep_pt = {};
    std::array<Float_t, ReadConst::MAX_GEN_LEP> gen_lep_eta = {};
    std::array<Float_t, ReadConst::MAX_GEN_LEP> gen_lep_phi = {};
    std::array<Float_t, ReadConst::MAX_GEN_LEP> gen_lep_mass = {};

    // 4x4x4 = 64
    std::array<Float_t, ReadConst::MAX_GEN_NU> gen_nu_pt = {};
    std::array<Float_t, ReadConst::MAX_GEN_NU> gen_nu_eta = {};
    std::array<Float_t, ReadConst::MAX_GEN_NU> gen_nu_phi = {};
    std::array<Float_t, ReadConst::MAX_GEN_NU> gen_nu_mass = {};

    // 4x4 = 16
    Float_t hbb_pt = 0.0f;
    Float_t hbb_eta = 0.0f;
    Float_t hbb_phi = 0.0f;
    Float_t hbb_mass = 0.0f;

    // 2x4 = 8
    Float_t reco_met_pt = 0.0f;
    Float_t reco_met_phi = 0.0f;

    // 8
    ULong64_t event_id = 0;

    // 4
    int n_reco_jet = 0;
    int n_reco_fatjet = 0;
}; 

#endif