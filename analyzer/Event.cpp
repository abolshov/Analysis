#include "Event.hpp"

Event::Event(TTree* tree, Channel ch) 
:   m_index(ch == Channel::SL ? GenTruthIdxMapSL : GenTruthIdxMapDL)
,   m_nu_index(ch == Channel::SL ? GenNuIdxMapSL : GenNuIdxMapDL)
,   m_reco_lep_index(ch == Channel::SL ? RecoLepIdxMapSL : RecoLepIdxMapDL)
,   m_branch_map(ch == Channel::SL ? GenTruthBranchMapSL : GenTruthBranchMapDL)
,   m_nu_branch_map(ch == Channel::SL ? GenNuBranchMapSL : GenNuBranchMapDL)
,   m_reco_lep_branch_map(ch == Channel::SL ? RecoLepBranchMapSL : RecoLepBranchMapDL)
,   gen_truth(m_index.size())
,   m_tree(tree)
{
    // // genjet data
    m_tree->SetBranchAddress("ncentralGenJet", &gen_jet.nGenJet);
    m_tree->SetBranchAddress("centralGenJet_pt", gen_jet.pt.get());
    m_tree->SetBranchAddress("centralGenJet_eta", gen_jet.eta.get());
    m_tree->SetBranchAddress("centralGenJet_phi", gen_jet.phi.get());
    m_tree->SetBranchAddress("centralGenJet_mass", gen_jet.mass.get());
    m_tree->SetBranchAddress("centralGenJet_partonFlavour", gen_jet.part_flav.get());
    m_tree->SetBranchAddress("centralGenJet_hadronFlavour", gen_jet.hadr_flav.get());

    // reco jet data
    m_tree->SetBranchAddress("ncentralJet", &reco_jet.nRecoJet);
    m_tree->SetBranchAddress("centralJet_pt", reco_jet.pt.get());
    m_tree->SetBranchAddress("centralJet_eta", reco_jet.eta.get());
    m_tree->SetBranchAddress("centralJet_phi", reco_jet.phi.get());
    m_tree->SetBranchAddress("centralJet_mass", reco_jet.mass.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", reco_jet.PNetRegPtRawCorr.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawRes", reco_jet.PNetRegPtRawRes.get());
    m_tree->SetBranchAddress("centralJet_btagPNetB", reco_jet.btagPNetB.get());
    m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", reco_jet.PNetRegPtRawCorrNu.get());

    // initialize gen truth variables
    for (auto const& [obj_name, branch_name]: m_branch_map)
    {
        size_t offset = m_index.at(obj_name);
        for (auto const& var_name: KinVarNames)
        {
            if (obj_name == "met" && (var_name == "_mass" || var_name == "_eta"))
            {
                continue;
            }

            std::string branch_full_name = branch_name + var_name;
            AddressFunc_t ptr = truth_address_map.at(var_name);
            Float_t* address = (this->*ptr)() + offset;
            m_tree->SetBranchAddress(branch_full_name.c_str(), address);
        }
    }

    // gen neutrino(s)
    for (auto const& [nu_name, branch_name]: m_nu_branch_map)
    {
        size_t offset = m_nu_index.at(nu_name);
        for (auto const& var_name: KinVarNames)
        {
            std::string branch_full_name = branch_name + var_name;
            AddressFunc_t ptr = nu_address_map.at(var_name);
            Float_t* address = (this->*ptr)() + offset;
            m_tree->SetBranchAddress(branch_full_name.c_str(), address);
        }

        m_tree->SetBranchAddress((branch_name + "_pdgId").c_str(), nu.pdgId.get() + offset);
    }

    // reco lepton(s)
    for (auto const& [lep_name, branch_name]: m_reco_lep_branch_map)
    {
        size_t offset = m_reco_lep_index.at(lep_name);
        for (auto const& var_name: KinVarNames)
        {
            std::string branch_full_name = branch_name + var_name;
            AddressFunc_t ptr = reco_lep_address_map.at(var_name);
            Float_t* address = (this->*ptr)() + offset;
            m_tree->SetBranchAddress(branch_full_name.c_str(), address);
        }

        m_tree->SetBranchAddress((branch_name + "_type").c_str(), reco_lep.lep_type.get() + offset);
        m_tree->SetBranchAddress((branch_name + "_iso").c_str(), reco_lep.lep_iso.get() + offset);
    }

    m_tree->SetBranchAddress("PuppiMET_pt", &reco_met_pt);
    m_tree->SetBranchAddress("PuppiMET_phi", &reco_met_phi);
}
