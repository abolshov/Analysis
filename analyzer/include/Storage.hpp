#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <array>

#include "TROOT.h"
#include "TTree.h"

#include "Constants.hpp"
#include "Definitions.hpp"

struct Storage
{   
    void ConnectTree(TTree* tree, Channel ch);

    //reco objects
    // 6x4x20 = 480
    std::array<Float_t, MAX_RECO_JET> reco_jet_pt = {};
    std::array<Float_t, MAX_RECO_JET> reco_jet_eta = {};
    std::array<Float_t, MAX_RECO_JET> reco_jet_phi = {};
    std::array<Float_t, MAX_RECO_JET> reco_jet_mass = {};
    std::array<Float_t, MAX_RECO_JET> reco_jet_corr = {};
    std::array<Float_t, MAX_RECO_JET> reco_jet_res = {};

    // 6x4x2 = 48
    std::array<Float_t, MAX_RECO_LEP> reco_lep_pt = {};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_eta = {};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_phi = {};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_mass = {};
    std::array<Int_t, MAX_RECO_LEP> reco_lep_type = {};
    std::array<Int_t, MAX_RECO_LEP> reco_lep_gen_kind = {};

    // 2x4 = 8
    Float_t reco_met_pt = 0.0;
    Float_t reco_met_phi = 0.0;

    // 8
    ULong64_t event_id = 0;

    // 4
    int n_reco_jet = 0;

    #ifdef DEBUG 
        // gen objects
        // 4x4x4 = 64
        std::array<Float_t, MAX_GEN_QUARK> gen_quark_pt = {};
        std::array<Float_t, MAX_GEN_QUARK> gen_quark_eta = {};
        std::array<Float_t, MAX_GEN_QUARK> gen_quark_phi = {};
        std::array<Float_t, MAX_GEN_QUARK> gen_quark_mass = {};

        // 4x4x4 = 64
        std::array<Float_t, MAX_GEN_LEP> gen_lep_pt = {};
        std::array<Float_t, MAX_GEN_LEP> gen_lep_eta = {};
        std::array<Float_t, MAX_GEN_LEP> gen_lep_phi = {};
        std::array<Float_t, MAX_GEN_LEP> gen_lep_mass = {};

        // 4x4x4 = 64
        std::array<Float_t, MAX_GEN_NU> gen_nu_pt = {};
        std::array<Float_t, MAX_GEN_NU> gen_nu_eta = {};
        std::array<Float_t, MAX_GEN_NU> gen_nu_phi = {};
        std::array<Float_t, MAX_GEN_NU> gen_nu_mass = {};

        // 4
        int n_gen_jet = 0;
    #endif 
}; 

#endif