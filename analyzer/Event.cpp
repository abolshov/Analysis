#include "Event.hpp"

Event::Event(TTree* tree) : m_tree(tree)
{
    m_tree->SetBranchAddress("ngenJet", &genjet.nGenJet);
    m_tree->SetBranchAddress("genJet_pt", genjet.pt.get());
    m_tree->SetBranchAddress("genJet_eta", genjet.eta.get());
    m_tree->SetBranchAddress("genJet_phi", genjet.phi.get());
    m_tree->SetBranchAddress("genJet_mass", genjet.mass.get());

    m_tree->SetBranchAddress("ncentralJet", &recojet.nRecoJet);
    m_tree->SetBranchAddress("centralJet_pt", recojet.pt.get());
    m_tree->SetBranchAddress("centralJet_eta", recojet.eta.get());
    m_tree->SetBranchAddress("centralJet_phi", recojet.phi.get());
    m_tree->SetBranchAddress("centralJet_mass", recojet.mass.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", recojet.PNetRegPtRawCorr.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawRes", recojet.PNetRegPtRawRes.get());
    m_tree->SetBranchAddress("centralJet_btagPNetB", recojet.btagPNetB.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", recojet.PNetRegPtRawCorrNu.get());

    m_tree->SetBranchAddress("PuppiMET_pt", &puppiMET_phi);
    m_tree->SetBranchAddress("PuppiMET_phi", &puppiMET_pt);
}