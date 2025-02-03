#include "SelectionUtils.hpp"
#include "MatchingTools.hpp"

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

bool CorrRecoLep(int lep_type, int lep_genLep_kind)
{
    bool reco_lep_mu = (lep_type == 2);
    bool reco_lep_ele = (lep_type == 1);

    bool gen_lep_mu = ((lep_genLep_kind == 2) || (lep_genLep_kind == 4));
    bool gen_lep_ele = ((lep_genLep_kind == 1) || (lep_genLep_kind == 3));

    bool corr_lep_reco = ((reco_lep_mu && gen_lep_mu) || (reco_lep_ele && gen_lep_ele));
    return corr_lep_reco;
}

bool IsRecoverable(Storage const& s, Channel ch, bool top_sel)
{
    // check all quarks
    int n_quarks = ch == Channel::DL ? 2 : 4;
    VecLVF_t quarks_p4;
    for (int i = 0; i < n_quarks; ++i)
    {
        Float_t max_eta = i < 2 ? 2.5 : 5.0;
        if (s.gen_quark_pt[i] < 20.0 || std::abs(s.gen_quark_eta[i]) >= max_eta)
        {
            return false;
        }
        quarks_p4.push_back(LorentzVectorF_t(s.gen_quark_pt[i], s.gen_quark_eta[i], s.gen_quark_phi[i], s.gen_quark_mass[i])); 
    }

    if (top_sel)
    {
        Float_t bb_dr = DeltaR(quarks_p4[static_cast<size_t>(Quark::b1)], quarks_p4[static_cast<size_t>(Quark::b2)]);
        if (bb_dr < 0.4)
        {
            return false;
        }

        if (ch == Channel::SL)
        {
            Float_t qq_dr = DeltaR(quarks_p4[static_cast<size_t>(Quark::q1)], quarks_p4[static_cast<size_t>(Quark::q2)]);
            if (qq_dr < 0.4)
            {
                return false;
            }
        }
    }

    // check all reco jets are present in the event
    int n_jets_min = ch == Channel::DL ? 2 : 4;
    int n_jets = s.n_reco_jet;
    if (n_jets < n_jets_min)
    {
        return false;
    }

    for (int i = 0; i < n_jets; ++i)
    {
        if (s.reco_jet_pt[i] < 20.0 || std::abs(s.reco_jet_eta[i]) >= 2.5)
        {
            return false;
        }
    }

    // check all reco leptons are present in the event
    int n_lep = ch == Channel::DL ? 2 : 1;
    for (int i = 0; i < n_lep; ++i)
    {
        if (s.reco_lep_pt[i] == 0.0)
        {
            return false;
        }
    }

    return true;
}

bool IsFiducial(Storage const& s, VecLVF_t const& jets, Channel ch)
{
    // check quark to reco jet matching
    int n_quarks = ch == Channel::DL ? 2 : 4;
    std::vector<LorentzVectorF_t> quarks_p4;
    for (int i = 0; i < n_quarks; ++i)
    {
        quarks_p4.push_back(LorentzVectorF_t(s.gen_quark_pt[i], s.gen_quark_eta[i], s.gen_quark_phi[i], s.gen_quark_mass[i]));
    }

    std::unordered_set<int> matches;
    for (auto const& qp4: quarks_p4)
    {
        int idx = MatchIdx(qp4, jets);
        matches.insert(idx);
    }

    if (matches.count(-1) || matches.size() != quarks_p4.size())
    {
        return false;
    }

    // check lepton reconstruction
    int n_lep = ch == Channel::DL ? 2 : 1;
    for (int i = 0; i < n_lep; ++i)
    {
        if (!CorrRecoLep(s.reco_lep_type[i], s.reco_lep_gen_kind[i]))
        {
            return false;
        }
    }

    return true;
}