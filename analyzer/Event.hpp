#ifndef EVENT_HPP
#define EVENT_HPP

#include "TTree.h"

#include "Constants.hpp"
#include "Objects.hpp"

class Event
{
    public:
    Event(TTree* tree);

    GenJet genjet;
    RecoJet recojet;
    
    Float_t puppiMET_phi;
    Float_t puppiMET_pt;

    private:
    TTree* m_tree;
};

#endif