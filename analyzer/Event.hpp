#ifndef EVENT_HPP
#define EVENT_HPP

#include "TTree.h"

#include "Constants.hpp"
#include "Objects.hpp"

class Event
{
    public:
    std::map<std::string, size_t> const m_index;
    std::map<std::string, size_t> const m_nu_index;
    std::map<std::string, size_t> const m_reco_lep_index;

    private:
    std::map<std::string, std::string> const m_branch_map;
    std::map<std::string, std::string> const m_nu_branch_map;
    std::map<std::string, std::string> const m_reco_lep_branch_map;

    public:
    Event(TTree* tree, Channel ch);

    GenJet_t gen_jet;        // gen level jets  
    RecoJet_t reco_jet;      // reco jets

    GenLep_t nu;             // gen level neutrinos
    RecoLep_t reco_lep;      // reco leptons

    Kinematics_t gen_truth;        // gen level quarks, lepton and met

    Float_t reco_met_pt;
    Float_t reco_met_phi;

    private:
    TTree* m_tree;

    typedef Float_t*(Event::*AddressFunc_t)();

    // truth addresses
    inline Float_t* GenTruthPt() { return gen_truth.pt.get(); }
    inline Float_t* GenTruthEta() { return gen_truth.eta.get(); }
    inline Float_t* GenTruthPhi() { return gen_truth.phi.get(); }
    inline Float_t* GenTruthMass() { return gen_truth.mass.get(); }

    inline static const std::map<std::string, AddressFunc_t> truth_address_map = { { "_pt", &Event::GenTruthPt },
                                                                                   { "_eta", &Event::GenTruthEta },
                                                                                   { "_phi", &Event::GenTruthPhi },
                                                                                   { "_mass", &Event::GenTruthMass } };

    // neutrino addresses
    inline Float_t* NuPt() { return nu.pt.get(); }
    inline Float_t* NuEta() { return nu.eta.get(); }
    inline Float_t* NuPhi() { return nu.phi.get(); }
    inline Float_t* NuMass() { return nu.mass.get(); }

    inline static const std::map<std::string, AddressFunc_t> nu_address_map = { { "_pt", &Event::NuPt },
                                                                                { "_eta", &Event::NuEta },
                                                                                { "_phi", &Event::NuPhi },
                                                                                { "_mass", &Event::NuMass } };

    // neutrino addresses
    inline Float_t* RecoLepPt() { return reco_lep.pt.get(); }
    inline Float_t* RecoLepEta() { return reco_lep.eta.get(); }
    inline Float_t* RecoLepPhi() { return reco_lep.phi.get(); }
    inline Float_t* RecoLepMass() { return reco_lep.mass.get(); }

    inline static const std::map<std::string, AddressFunc_t> reco_lep_address_map = { { "_pt", &Event::RecoLepPt },
                                                                                      { "_eta", &Event::RecoLepEta },
                                                                                      { "_phi", &Event::RecoLepPhi },
                                                                                      { "_mass", &Event::RecoLepMass } };
};

#endif