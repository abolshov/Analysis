#ifndef EVENT_HPP
#define EVENT_HPP

#include "TTree.h"

#include "Constants.hpp"
#include "Objects.hpp"
#include "Input.hpp"

using GenJet_t = GenJet;
using RecoJet_t = RecoJet;
using Particle_t = Particle;
using Kin_t = Kinematics;
using GenLep_t = GenLep;
using RecoLep_t = RecoLep;

class Event
{
    public:
    Event(TTree* tree);

    GenJet_t genjet;        // gen level jets  
    RecoJet_t recojet;      // reco jets

    Kin_t gen_truth;        // gen level quarks, lepton and met
    GenLep_t gen_lep;

    Particle_t nu;
    
    Particle_t reco_met;
    Particle_t reco_lep;

    EstimatorInput MakeEstimatorInput(std::string const& pdf_file_name) const;
    ValidatorInput MakeValidatorInput(std::string const& pdf_file_name) const;

    private:
    TTree* m_tree;

    inline Float_t* GenTruthPt() { return gen_truth.pt.get(); }
    inline Float_t* GenTruthEta() { return gen_truth.eta.get(); }
    inline Float_t* GenTruthPhi() { return gen_truth.phi.get(); }
    inline Float_t* GenTruthMass() { return gen_truth.mass.get(); }

    typedef Float_t*(Event::*AddressFunc_t)();
    inline static const std::map<std::string, AddressFunc_t> address_map = { { "_pt", &Event::GenTruthPt },
                                                                             { "_eta", &Event::GenTruthEta },
                                                                             { "_phi", &Event::GenTruthPhi },
                                                                             { "_mass", &Event::GenTruthMass } };
    

    inline size_t Offset(Obj obj) { return static_cast<size_t>(obj); }
};

#endif