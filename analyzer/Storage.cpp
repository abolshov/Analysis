#include "Storage.hpp"

Storage::Storage(TTree* tree, Channel ch)
:   m_tree(tree)
{   
    // reco jet data
    // m_tree->SetBranchAddress("ncentralJet", &reco_jet.nRecoJet);
    m_tree->SetBranchAddress("centralJet_pt", reco_jet_pt.data());
    m_tree->SetBranchAddress("centralJet_eta", reco_jet_eta.data());
    m_tree->SetBranchAddress("centralJet_phi", reco_jet_phi.data());
    m_tree->SetBranchAddress("centralJet_mass", reco_jet_mass.data());
    // m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", reco_jet.PNetRegPtRawCorr.get());
    // m_tree->SetBranchAddress("centralJet_PNetRegPtRawRes", reco_jet.PNetRegPtRawRes.get());
    // m_tree->SetBranchAddress("centralJet_btagPNetB", reco_jet.btagPNetB.get());
    // m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", reco_jet.PNetRegPtRawCorrNu.get());

    m_tree->SetBranchAddress("lep1_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("lep1_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("lep1_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("lep1_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep1));

    m_tree->SetBranchAddress("lep2_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep2));
    m_tree->SetBranchAddress("lep2_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep2));
    m_tree->SetBranchAddress("lep2_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep2));
    m_tree->SetBranchAddress("lep2_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep2));

    // genjet data
    // m_tree->SetBranchAddress("ncentralGenJet", &gen_jet.nGenJet);
    m_tree->SetBranchAddress("centralGenJet_pt", gen_jet_pt.data());
    m_tree->SetBranchAddress("centralGenJet_eta", gen_jet_eta.data());
    m_tree->SetBranchAddress("centralGenJet_phi", gen_jet_phi.data());
    m_tree->SetBranchAddress("centralGenJet_mass", gen_jet_mass.data());

    m_tree->SetBranchAddress("genV1prod1_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("genV1prod1_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("genV1prod1_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep1));
    m_tree->SetBranchAddress("genV1prod1_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep1));

    if (ch == Channel::DL)
    {
        m_tree->SetBranchAddress("genV2prod1_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep2));
        m_tree->SetBranchAddress("genV2prod1_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep2));
        m_tree->SetBranchAddress("genV2prod1_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep2));
        m_tree->SetBranchAddress("genV2prod1_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep2));
    }

    m_tree->SetBranchAddress("genb1_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::b1));
    m_tree->SetBranchAddress("genb1_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::b1));
    m_tree->SetBranchAddress("genb1_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::b1));
    m_tree->SetBranchAddress("genb1_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::b1));

    m_tree->SetBranchAddress("genb2_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::b2));
    m_tree->SetBranchAddress("genb2_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::b2));
    m_tree->SetBranchAddress("genb2_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::b2));
    m_tree->SetBranchAddress("genb2_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::b2));

    if (ch == Channel::SL)
    {
        m_tree->SetBranchAddress("genV2prod1_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::q1));
        m_tree->SetBranchAddress("genV2prod1_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::q1));
        m_tree->SetBranchAddress("genV2prod1_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::q1));
        m_tree->SetBranchAddress("genV2prod1_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::q1));

        m_tree->SetBranchAddress("genV2prod2_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::q2));
        m_tree->SetBranchAddress("genV2prod2_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::q2));
        m_tree->SetBranchAddress("genV2prod2_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::q2));
        m_tree->SetBranchAddress("genV2prod2_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::q2));
    }
}