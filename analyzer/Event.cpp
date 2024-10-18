#include "Event.hpp"

Event::Event(TTree* tree, Channel ch) 
: m_index(ch == Channel::SL ? GenTruthIdxMapSL : GenTruthIdxMapDL),
  m_branch_map(ch == Channel::SL ? GenTruthBranchMapSL : GenTruthBranchMapDL),
  gen_truth(m_index.size()),
  m_tree(tree)
{
    // // genjet data
    m_tree->SetBranchAddress("ncentralGenJet", &genjet.nGenJet);
    m_tree->SetBranchAddress("centralGenJet_pt", genjet.pt.get());
    m_tree->SetBranchAddress("centralGenJet_eta", genjet.eta.get());
    m_tree->SetBranchAddress("centralGenJet_phi", genjet.phi.get());
    m_tree->SetBranchAddress("centralGenJet_mass", genjet.mass.get());
    m_tree->SetBranchAddress("centralGenJet_partonFlavour", genjet.part_flav.get());
    m_tree->SetBranchAddress("centralGenJet_hadronFlavour", genjet.hadr_flav.get());

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

    // initialize gen truth variables
    for (auto const& [obj_name, branch_name]: m_branch_map)
    {
        for (auto const& var_name: KinVarNames)
        {
            if (obj_name == "met" && (var_name == "_mass" || var_name == "_eta"))
            {
                continue;
            }

            std::string branch_full_name = branch_name + var_name;
            size_t offset = m_index.at(obj_name);
            AddressFunc_t ptr = address_map.at(var_name);
            Float_t* address = (this->*ptr)() + offset;

            std::cout << obj_name << ", " << branch_full_name << ", " << address << "\n";

            m_tree->SetBranchAddress(branch_full_name.c_str(), address);
        }
    }
}

// Event::Event(TTree* tree) : gen_truth(static_cast<size_t>(Obj::count)), m_tree(tree)
// {
//     // genjet data
//     m_tree->SetBranchAddress("ncentralGenJet", &genjet.nGenJet);
//     m_tree->SetBranchAddress("centralGenJet_pt", genjet.pt.get());
//     m_tree->SetBranchAddress("centralGenJet_eta", genjet.eta.get());
//     m_tree->SetBranchAddress("centralGenJet_phi", genjet.phi.get());
//     m_tree->SetBranchAddress("centralGenJet_mass", genjet.mass.get());
//     m_tree->SetBranchAddress("centralGenJet_partonFlavour", genjet.part_flav.get());
//     m_tree->SetBranchAddress("centralGenJet_hadronFlavour", genjet.hadr_flav.get());

//     // true gen partons data
//     m_tree->SetBranchAddress("genb1_pt", gen_truth.pt.get() + Offset(Obj::b1));
//     m_tree->SetBranchAddress("genb2_pt", gen_truth.pt.get() + Offset(Obj::b2));
//     m_tree->SetBranchAddress("genV2prod1_pt", gen_truth.pt.get() + Offset(Obj::q1));
//     m_tree->SetBranchAddress("genV2prod2_pt", gen_truth.pt.get() + Offset(Obj::q2));
//     m_tree->SetBranchAddress("genV1prod1_pt", gen_truth.pt.get() + Offset(Obj::lep));

//     m_tree->SetBranchAddress("genb1_eta", gen_truth.eta.get() + Offset(Obj::b1));
//     m_tree->SetBranchAddress("genb2_eta", gen_truth.eta.get() + Offset(Obj::b2));
//     m_tree->SetBranchAddress("genV2prod1_eta", gen_truth.eta.get() + Offset(Obj::q1));
//     m_tree->SetBranchAddress("genV2prod2_eta", gen_truth.eta.get() + Offset(Obj::q2));
//     m_tree->SetBranchAddress("genV1prod1_eta", gen_truth.eta.get() + Offset(Obj::lep));

//     m_tree->SetBranchAddress("genb1_phi", gen_truth.phi.get() + Offset(Obj::b1));
//     m_tree->SetBranchAddress("genb2_phi", gen_truth.phi.get() + Offset(Obj::b2));
//     m_tree->SetBranchAddress("genV2prod1_phi", gen_truth.phi.get() + Offset(Obj::q1));
//     m_tree->SetBranchAddress("genV2prod2_phi", gen_truth.phi.get() + Offset(Obj::q2));
//     m_tree->SetBranchAddress("genV1prod1_phi", gen_truth.phi.get() + Offset(Obj::lep));

//     m_tree->SetBranchAddress("genb1_mass", gen_truth.mass.get() + Offset(Obj::b1));
//     m_tree->SetBranchAddress("genb2_mass", gen_truth.mass.get() + Offset(Obj::b2));
//     m_tree->SetBranchAddress("genV2prod1_mass", gen_truth.mass.get() + Offset(Obj::q1));
//     m_tree->SetBranchAddress("genV2prod2_mass", gen_truth.mass.get() + Offset(Obj::q2));
//     m_tree->SetBranchAddress("genV1prod1_mass", gen_truth.mass.get() + Offset(Obj::lep));

//     m_tree->SetBranchAddress("GenMET_pt", gen_truth.pt.get() + Offset(Obj::met));
//     m_tree->SetBranchAddress("GenMET_phi", gen_truth.phi.get() + Offset(Obj::met));

//     // true neutrino
//     m_tree->SetBranchAddress("genV1prod2_pt", &nu.pt);
//     m_tree->SetBranchAddress("genV1prod2_eta", &nu.eta);
//     m_tree->SetBranchAddress("genV1prod2_phi", &nu.phi);
//     m_tree->SetBranchAddress("genV1prod2_phi", &nu.mass);

//     // reco jet data
//     m_tree->SetBranchAddress("ncentralJet", &recojet.nRecoJet);
//     m_tree->SetBranchAddress("centralJet_pt", recojet.pt.get());
//     m_tree->SetBranchAddress("centralJet_eta", recojet.eta.get());
//     m_tree->SetBranchAddress("centralJet_phi", recojet.phi.get());
//     m_tree->SetBranchAddress("centralJet_mass", recojet.mass.get());
//     m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorr", recojet.PNetRegPtRawCorr.get());
//     m_tree->SetBranchAddress("centralJet_PNetRegPtRawRes", recojet.PNetRegPtRawRes.get());
//     m_tree->SetBranchAddress("centralJet_btagPNetB", recojet.btagPNetB.get());
//     m_tree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", recojet.PNetRegPtRawCorrNu.get());

//     m_tree->SetBranchAddress("PuppiMET_pt", &reco_met.pt);
//     m_tree->SetBranchAddress("PuppiMET_phi", &reco_met.phi);

//     m_tree->SetBranchAddress("lep1_pt", &reco_lep.pt);
//     m_tree->SetBranchAddress("lep1_eta", &reco_lep.eta);
//     m_tree->SetBranchAddress("lep1_phi", &reco_lep.phi);
//     m_tree->SetBranchAddress("lep1_mass", &reco_lep.mass);
// }

EstimatorInput Event::MakeEstimatorInput(std::string const& pdf_file_name) const
{
    std::vector<TLorentzVector> p4;
    size_t n_inputs = static_cast<size_t>(Obj::count);
    p4.reserve(n_inputs);

    return EstimatorInput(std::move(p4), pdf_file_name);
}