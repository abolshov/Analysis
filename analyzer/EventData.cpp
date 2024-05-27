#include "EventData.hpp"

EventData::EventData(TTree& tree) : m_tree(tree) 
{
    m_tree.SetBranchAddress("nGenPart", &nGenPart);
    m_tree.SetBranchAddress("GenPart_eta", &GenPart_eta);
    m_tree.SetBranchAddress("GenPart_mass", &GenPart_mass);
    m_tree.SetBranchAddress("GenPart_phi", &GenPart_phi);
    m_tree.SetBranchAddress("GenPart_pt", &GenPart_pt);
    m_tree.SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    m_tree.SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    m_tree.SetBranchAddress("GenPart_status", &GenPart_status);
    m_tree.SetBranchAddress("GenPart_phi", &GenPart_phi);

    m_tree.SetBranchAddress("nGenJet", &nGenJetAK4);
    m_tree.SetBranchAddress("GenJet_eta", &GenJetAK4_eta);
    m_tree.SetBranchAddress("GenJet_mass", &GenJetAK4_mass);
    m_tree.SetBranchAddress("GenJet_phi", &GenJetAK4_phi);
    m_tree.SetBranchAddress("GenJet_pt", &GenJetAK4_pt);
    m_tree.SetBranchAddress("GenJet_partonFlavour", &GenJetAK4_partonFlavour);
    m_tree.SetBranchAddress("GenJet_hadronFlavour", &GenJetAK4_hadronFlavour);

    m_tree.SetBranchAddress("GenMET_phi", &GenMET_phi);
    m_tree.SetBranchAddress("GenMET_pt", &GenMET_pt);
}