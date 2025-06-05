#include "Buffer.hpp"

void Buffer::ConnectTree(TTree* tree, Channel ch)
{   
    // reco jet data
    tree->SetBranchAddress("ncentralJet", &n_reco_jet);
    tree->SetBranchAddress("centralJet_pt", reco_jet_pt.data());
    tree->SetBranchAddress("centralJet_eta", reco_jet_eta.data());
    tree->SetBranchAddress("centralJet_phi", reco_jet_phi.data());
    tree->SetBranchAddress("centralJet_mass", reco_jet_mass.data());
    tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", reco_jet_corr.data());
    tree->SetBranchAddress("centralJet_PNetRegPtRawRes", reco_jet_res.data());
    tree->SetBranchAddress("centralJet_btagPNetB", reco_jet_btag.data());
    tree->SetBranchAddress("centralJet_btagPNetQvG", reco_jet_qvg.data());

    // reco fatjet data
    tree->SetBranchAddress("ncentralJet", &n_reco_fatjet);
    tree->SetBranchAddress("SelectedFatJet_pt", reco_fatjet_pt.data());
    tree->SetBranchAddress("SelectedFatJet_eta", reco_fatjet_eta.data());
    tree->SetBranchAddress("SelectedFatJet_phi", reco_fatjet_phi.data());
    tree->SetBranchAddress("SelectedFatJet_mass", reco_fatjet_mass.data());
    tree->SetBranchAddress("SelectedFatJet_particleNet_HbbvsQCD", reco_fatjet_btag_HbbvsQCD.data());

    tree->SetBranchAddress("lep1_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("lep1_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("lep1_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("lep1_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("lep1_type", reco_lep_type.data() + static_cast<size_t>(Lep::lep1));

    tree->SetBranchAddress("lep2_pt", reco_lep_pt.data() + static_cast<size_t>(Lep::lep2));
    tree->SetBranchAddress("lep2_eta", reco_lep_eta.data() + static_cast<size_t>(Lep::lep2));
    tree->SetBranchAddress("lep2_phi", reco_lep_phi.data() + static_cast<size_t>(Lep::lep2));
    tree->SetBranchAddress("lep2_mass", reco_lep_mass.data() + static_cast<size_t>(Lep::lep2));
    tree->SetBranchAddress("lep2_type", reco_lep_type.data() + static_cast<size_t>(Lep::lep2));

    #ifdef OLD_LEP_BRANCH
        tree->SetBranchAddress("lep1_genLep_kind", reco_lep_gen_kind.data() + static_cast<size_t>(Lep::lep1));
        tree->SetBranchAddress("lep2_genLep_kind", reco_lep_gen_kind.data() + static_cast<size_t>(Lep::lep2));
    #else 
        tree->SetBranchAddress("lep1_gen_kind", reco_lep_gen_kind.data() + static_cast<size_t>(Lep::lep1));
        tree->SetBranchAddress("lep2_gen_kind", reco_lep_gen_kind.data() + static_cast<size_t>(Lep::lep2));
    #endif

    tree->SetBranchAddress("event", &event_id);
    tree->SetBranchAddress("PuppiMET_pt", &reco_met_pt);
    tree->SetBranchAddress("PuppiMET_phi", &reco_met_phi);

    tree->SetBranchAddress("genV1prod1_pt", gen_lep_pt.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("genV1prod1_eta", gen_lep_eta.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("genV1prod1_phi", gen_lep_phi.data() + static_cast<size_t>(Lep::lep1));
    tree->SetBranchAddress("genV1prod1_mass", gen_lep_mass.data() + static_cast<size_t>(Lep::lep1));

    if (ch == Channel::DL)
    {
        tree->SetBranchAddress("genV2prod1_pt", gen_lep_pt.data() + static_cast<size_t>(Lep::lep2));
        tree->SetBranchAddress("genV2prod1_eta", gen_lep_eta.data() + static_cast<size_t>(Lep::lep2));
        tree->SetBranchAddress("genV2prod1_phi", gen_lep_phi.data() + static_cast<size_t>(Lep::lep2));
        tree->SetBranchAddress("genV2prod1_mass", gen_lep_mass.data() + static_cast<size_t>(Lep::lep2));
    }

    tree->SetBranchAddress("genb1_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::b1));
    tree->SetBranchAddress("genb1_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::b1));
    tree->SetBranchAddress("genb1_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::b1));
    tree->SetBranchAddress("genb1_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::b1));

    tree->SetBranchAddress("genb2_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::b2));
    tree->SetBranchAddress("genb2_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::b2));
    tree->SetBranchAddress("genb2_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::b2));
    tree->SetBranchAddress("genb2_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::b2));

    if (ch == Channel::SL)
    {
        tree->SetBranchAddress("genV2prod1_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::q1));
        tree->SetBranchAddress("genV2prod1_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::q1));
        tree->SetBranchAddress("genV2prod1_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::q1));
        tree->SetBranchAddress("genV2prod1_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::q1));

        tree->SetBranchAddress("genV2prod2_pt", gen_quark_pt.data() + static_cast<size_t>(Quark::q2));
        tree->SetBranchAddress("genV2prod2_eta", gen_quark_eta.data() + static_cast<size_t>(Quark::q2));
        tree->SetBranchAddress("genV2prod2_phi", gen_quark_phi.data() + static_cast<size_t>(Quark::q2));
        tree->SetBranchAddress("genV2prod2_mass", gen_quark_mass.data() + static_cast<size_t>(Quark::q2));
    }

    tree->SetBranchAddress("genHbb_pt", &hbb_pt);
    tree->SetBranchAddress("genHbb_eta", &hbb_eta);
    tree->SetBranchAddress("genHbb_phi", &hbb_phi);
    tree->SetBranchAddress("genHbb_mass", &hbb_mass);
}