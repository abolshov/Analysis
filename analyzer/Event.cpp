#include "Event.hpp"

Event::Event(TTree* tree) : gen_truth(ObjSLRes::count), m_tree(tree)
{
    // genjet data
    m_tree->SetBranchAddress("ncentralGenJet", &genjet.nGenJet);
    m_tree->SetBranchAddress("centralGenJet_pt", genjet.pt.get());
    m_tree->SetBranchAddress("centralGenJet_eta", genjet.eta.get());
    m_tree->SetBranchAddress("centralGenJet_phi", genjet.phi.get());
    m_tree->SetBranchAddress("centralGenJet_mass", genjet.mass.get());
    m_tree->SetBranchAddress("centralGenJet_partonFlavour", genjet.part_flav.get());
    m_tree->SetBranchAddress("centralGenJet_hadronFlavour", genjet.hadr_flav.get());

    // true gen partons data
    m_tree->SetBranchAddress("genb1_pt", gen_truth.pt.get() + ObjSLRes::b1);
    m_tree->SetBranchAddress("genb2_pt", gen_truth.pt.get() + ObjSLRes::b2);
    m_tree->SetBranchAddress("genV2prod1_pt", gen_truth.pt.get() + ObjSLRes::q1);
    m_tree->SetBranchAddress("genV2prod2_pt", gen_truth.pt.get() + ObjSLRes::q2);
    m_tree->SetBranchAddress("genV1prod1_pt", gen_truth.pt.get() + ObjSLRes::lep);

    m_tree->SetBranchAddress("genb1_eta", gen_truth.eta.get() + ObjSLRes::b1);
    m_tree->SetBranchAddress("genb2_eta", gen_truth.eta.get() + ObjSLRes::b2);
    m_tree->SetBranchAddress("genV2prod1_eta", gen_truth.eta.get() + ObjSLRes::q1);
    m_tree->SetBranchAddress("genV2prod2_eta", gen_truth.eta.get() + ObjSLRes::q2);
    m_tree->SetBranchAddress("genV1prod1_eta", gen_truth.eta.get() + ObjSLRes::lep);

    m_tree->SetBranchAddress("genb1_phi", gen_truth.phi.get() + ObjSLRes::b1);
    m_tree->SetBranchAddress("genb2_phi", gen_truth.phi.get() + ObjSLRes::b2);
    m_tree->SetBranchAddress("genV2prod1_phi", gen_truth.phi.get() + ObjSLRes::q1);
    m_tree->SetBranchAddress("genV2prod2_phi", gen_truth.phi.get() + ObjSLRes::q2);
    m_tree->SetBranchAddress("genV1prod1_phi", gen_truth.phi.get() + ObjSLRes::lep);

    m_tree->SetBranchAddress("genb1_mass", gen_truth.mass.get() + ObjSLRes::b1);
    m_tree->SetBranchAddress("genb2_mass", gen_truth.mass.get() + ObjSLRes::b2);
    m_tree->SetBranchAddress("genV2prod1_mass", gen_truth.mass.get() + ObjSLRes::q1);
    m_tree->SetBranchAddress("genV2prod2_mass", gen_truth.mass.get() + ObjSLRes::q2);
    m_tree->SetBranchAddress("genV1prod1_mass", gen_truth.mass.get() + ObjSLRes::lep);

    m_tree->SetBranchAddress("GenMET_pt", gen_truth.pt.get() + ObjSLRes::met);
    m_tree->SetBranchAddress("GenMET_phi", gen_truth.phi.get() + ObjSLRes::met);

    // true neutrino
    m_tree->SetBranchAddress("genV1prod2_pt", &nu.pt);
    m_tree->SetBranchAddress("genV1prod2_eta", &nu.eta);
    m_tree->SetBranchAddress("genV1prod2_phi", &nu.phi);
    m_tree->SetBranchAddress("genV1prod2_phi", &nu.mass);

    // reco jet data
    m_tree->SetBranchAddress("ncentralJet", &recojet.nRecoJet);
    m_tree->SetBranchAddress("centralJet_pt", recojet.pt.get());
    m_tree->SetBranchAddress("centralJet_eta", recojet.eta.get());
    m_tree->SetBranchAddress("centralJet_phi", recojet.phi.get());
    m_tree->SetBranchAddress("centralJet_mass", recojet.mass.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", recojet.PNetRegPtRawCorr.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawRes", recojet.PNetRegPtRawRes.get());
    m_tree->SetBranchAddress("centralJet_btagPNetB", recojet.btagPNetB.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", recojet.PNetRegPtRawCorrNu.get());

    m_tree->SetBranchAddress("PuppiMET_pt", &reco_met.pt);
    m_tree->SetBranchAddress("PuppiMET_phi", &reco_met.phi);

    m_tree->SetBranchAddress("lep1_pt", &reco_lep.pt);
    m_tree->SetBranchAddress("lep1_eta", &reco_lep.eta);
    m_tree->SetBranchAddress("lep1_phi", &reco_lep.phi);
    m_tree->SetBranchAddress("lep1_mass", &reco_lep.mass);
}