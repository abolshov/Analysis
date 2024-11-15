#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <array>

#include "TROOT.h"

#include "Constants.hpp"

struct Storage
{   
    Storage(TTree* tree, Channel ch);
    //reco objects
    std::array<Float_t, MAX_RECO_JET> reco_jet_pt = {0.0};
    std::array<Float_t, MAX_RECO_JET> reco_jet_eta = {0.0};
    std::array<Float_t, MAX_RECO_JET> reco_jet_phi = {0.0};
    std::array<Float_t, MAX_RECO_JET> reco_jet_mass = {0.0};

    std::array<Float_t, MAX_RECO_LEP> reco_lep_pt = {0.0};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_eta = {0.0};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_phi = {0.0};
    std::array<Float_t, MAX_RECO_LEP> reco_lep_mass = {0.0};

    // gen objects
    std::array<Float_t, MAX_GEN_JET> gen_jet_pt = {0.0};
    std::array<Float_t, MAX_GEN_JET> gen_jet_eta = {0.0};
    std::array<Float_t, MAX_GEN_JET> gen_jet_phi = {0.0};
    std::array<Float_t, MAX_GEN_JET> gen_jet_mass = {0.0};

    std::array<Float_t, MAX_GEN_LEP> gen_lep_pt = {0.0};
    std::array<Float_t, MAX_GEN_LEP> gen_lep_eta = {0.0};
    std::array<Float_t, MAX_GEN_LEP> gen_lep_phi = {0.0};
    std::array<Float_t, MAX_GEN_LEP> gen_lep_mass = {0.0};

    std::array<Float_t, MAX_GEN_NU> gen_nu_pt = {0.0};
    std::array<Float_t, MAX_GEN_NU> gen_nu_eta = {0.0};
    std::array<Float_t, MAX_GEN_NU> gen_nu_phi = {0.0};
    std::array<Float_t, MAX_GEN_NU> gen_nu_mass = {0.0};

    std::array<Float_t, MAX_GEN_QUARK> gen_quark_pt = {0.0};
    std::array<Float_t, MAX_GEN_QUARK> gen_quark_eta = {0.0};
    std::array<Float_t, MAX_GEN_QUARK> gen_quark_phi = {0.0};
    std::array<Float_t, MAX_GEN_QUARK> gen_quark_mass = {0.0};

    private:
    TTree* m_tree;
}; 

#endif