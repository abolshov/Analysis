#include "EstimatorDoubleLep.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_set>

#include "TVector2.h"
#include "Math/GenVector/VectorUtil.h" // DeltaPhi
using ROOT::Math::VectorUtil::DeltaPhi;

#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"

EstimatorDoubleLep::EstimatorDoubleLep(TString const& pdf_file_name, AggregationMode aggr_mode)
:   EstimatorBase(aggr_mode)
,   m_iter_data(std::make_unique<IterData>())
{
    m_pdf_1d.resize(pdf1d_dl_names.size());
    m_pdf_2d.resize(pdf2d_dl_names.size());
    
    TFile* pf = TFile::Open(pdf_file_name);
    Get1dPDFs(pf, m_pdf_1d, Channel::DL);
    Get2dPDFs(pf, m_pdf_2d, Channel::DL);
    pf->Close();
}

ArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateCombination(VecLVF_t const& particles, ULong64_t evt_id, TString const& comb_label)
{
    ArrF_t<ESTIM_OUT_SZ> res{};
    std::fill(res.begin(), res.end(), -1.0f);

    LorentzVectorF_t const& bj1 = particles[static_cast<size_t>(ObjDL::bj1)];
    LorentzVectorF_t const& bj2 = particles[static_cast<size_t>(ObjDL::bj2)];
    LorentzVectorF_t const& lep1 = particles[static_cast<size_t>(ObjDL::lep1)];
    LorentzVectorF_t const& lep2 = particles[static_cast<size_t>(ObjDL::lep2)];
    LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjDL::met)];

    UHist_t<TH1F>& pdf_b1 = m_pdf_1d[static_cast<size_t>(PDF1_dl::b1)];
    UHist_t<TH1F>& pdf_mw_onshell = m_pdf_1d[static_cast<size_t>(PDF1_dl::mw_onshell)];

    if (m_aggr_mode == AggregationMode::Combination)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt_id, comb_label.Data()));
    }

    [[maybe_unused]] TString tree_name = Form("evt_%llu_%s", evt_id, comb_label.Data());
    if (m_recorder.ShouldRecord())
    {
        m_recorder.ResetTree(MakeTree(tree_name));
    }

    for (int i = 0; i < N_ITER; ++i)
    {
        Float_t eta_gen = m_prg->Uniform(-6, 6);
        Float_t phi_gen = m_prg->Uniform(-3.1415926, 3.1415926);
        Float_t mh = m_prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);
        Float_t mw = pdf_mw_onshell->GetRandom(m_prg.get());
        Float_t smear_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t smear_dpy = m_prg->Gaus(0.0, MET_SIGMA);

        auto bresc = ComputeJetResc(bj1, bj2, pdf_b1, mh);
        if (!bresc.has_value())
        {
            continue;
        }
        auto [c1, c2] = bresc.value();

        LorentzVectorF_t b1 = bj1;
        LorentzVectorF_t b2 = bj2;
        b1 *= c1;
        b2 *= c2;

        Float_t bjet_resc_dpx = -1.0*(c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px();
        Float_t bjet_resc_dpy = -1.0*(c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py();

        Float_t met_corr_px = met.Px() + bjet_resc_dpx + smear_dpx;
        Float_t met_corr_py = met.Py() + bjet_resc_dpy + smear_dpy;

        Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        LorentzVectorF_t met_corr(met_corr_pt, 0.0, met_corr_phi, 0.0);

        LorentzVectorF_t Hbb = b1;
        Hbb += b2;

        std::vector<Float_t> estimates;
        std::array<LorentzVectorF_t, CONTROL> onshellW{};
        std::array<LorentzVectorF_t, CONTROL> offshellW{};
        std::array<LorentzVectorF_t, CONTROL> l_offshell{};
        std::array<LorentzVectorF_t, CONTROL> l_onshell{};
        std::array<LorentzVectorF_t, CONTROL> Hww{};
        std::array<LorentzVectorF_t, CONTROL> Xhh{};
        std::array<LorentzVectorF_t, CONTROL> nu_offshell{};
        std::array<LorentzVectorF_t, CONTROL> nu_onshell{};
        std::array<Float_t, CONTROL> mass{};
        // two options: 
        // lep1 comes from onshell W and lep2 comes from offshell W
        // lep1 comes from offshell W and lep2 comes from onshell W
        // when neutrino is computed in each case again two options are possible: delta_eta is added or subtracted to nu
        // in total: 4 combinations; they are encoded in this for loop
        for (int control = 0; control < CONTROL; ++control)
        {
            int is_onshell = control / 2;
            if (is_onshell == 0) 
            {
                l_onshell[control] = lep1;
                l_offshell[control] = lep2;
            } 
            else 
            {
                l_onshell[control] = lep2;
                l_offshell[control] = lep1;
            }

            auto opt_nu_onshell = NuFromOnshellW(eta_gen, phi_gen, mw, l_onshell[control]);  
            if (!opt_nu_onshell)
            {
                continue;
            }
            nu_onshell[control] = opt_nu_onshell.value();

            int is_offshell = control % 2;
            auto opt_nu_offshell = NuFromOffshellW(lep1, lep2, nu_onshell[control], met_corr, is_offshell, mh);
            if (!opt_nu_offshell)
            {
                continue;
            }
            nu_offshell[control] = opt_nu_offshell.value();

            onshellW[control] = l_onshell[control] + nu_onshell[control];
            offshellW[control] = l_offshell[control] + nu_offshell[control];
            Hww[control] = onshellW[control] + offshellW[control];

            if (offshellW[control].M() > mh/2)
            {
                continue;
            }

            if (std::abs(Hww[control].M() - mh) > 1.0)
            {
                continue;
            }

            Xhh[control] = Hbb + Hww[control];
            mass[control] = Xhh[control].M();
            if (mass[control] > 0.0)
            {
                estimates.push_back(mass[control]);
            }
        }

        if (estimates.empty())
        {
            continue;
        }

        Int_t num_sol = estimates.size();
        Float_t weight = estimates.empty() ? 0.0 : 1.0/estimates.size();
        for (auto est: estimates)
        {
            m_res_mass->Fill(est, weight);
        }

        if (m_recorder.ShouldRecord())
        {
            for (int control = 0; control < CONTROL; ++control)
            {
                m_iter_data->offshellW_pt[control] = offshellW[control].Pt();
                m_iter_data->offshellW_eta[control] = offshellW[control].Eta();
                m_iter_data->offshellW_phi[control] = offshellW[control].Phi();
                m_iter_data->offshellW_mass[control] = offshellW[control].M();

                m_iter_data->onshellW_pt[control] = onshellW[control].Pt();
                m_iter_data->onshellW_eta[control] = onshellW[control].Eta();
                m_iter_data->onshellW_phi[control] = onshellW[control].Phi();
                m_iter_data->onshellW_mass[control] = onshellW[control].M();

                m_iter_data->l_offshell_pt[control] = l_offshell[control].Pt();
                m_iter_data->l_offshell_eta[control] = l_offshell[control].Eta();
                m_iter_data->l_offshell_phi[control] = l_offshell[control].Phi();
                m_iter_data->l_offshell_mass[control] = l_offshell[control].M();

                m_iter_data->l_onshell_pt[control] = l_onshell[control].Pt();
                m_iter_data->l_onshell_eta[control] = l_onshell[control].Eta();
                m_iter_data->l_onshell_phi[control] = l_onshell[control].Phi();
                m_iter_data->l_onshell_mass[control] = l_onshell[control].M();

                m_iter_data->Hww_pt[control] = Hww[control].Pt();
                m_iter_data->Hww_eta[control] = Hww[control].Eta();
                m_iter_data->Hww_phi[control] = Hww[control].Phi();
                m_iter_data->Hww_mass[control] = Hww[control].M();

                m_iter_data->nu_offshell_pt[control] = nu_offshell[control].Pt();
                m_iter_data->nu_offshell_eta[control] = nu_offshell[control].Eta();
                m_iter_data->nu_offshell_phi[control] = nu_offshell[control].Phi();

                m_iter_data->nu_onshell_pt[control] = nu_onshell[control].Pt();
                m_iter_data->nu_onshell_eta[control] = nu_onshell[control].Eta();
                m_iter_data->nu_onshell_phi[control] = nu_onshell[control].Phi();

                m_iter_data->Xhh_pt[control] = Xhh[control].Pt();
                m_iter_data->Xhh_eta[control] = Xhh[control].Eta();
                m_iter_data->Xhh_phi[control] = Xhh[control].Phi();
                m_iter_data->Xhh_mass[control] = Xhh[control].M();

                m_iter_data->mass[control] = mass[control];
            }

            m_iter_data->met_corr_pt = met_corr.Pt();
            m_iter_data->met_corr_phi = met_corr.Phi();

            m_iter_data->b1_pt = b1.Pt();
            m_iter_data->b1_eta = b1.Eta();
            m_iter_data->b1_phi = b1.Phi();
            m_iter_data->b1_mass = b1.M();

            m_iter_data->b2_pt = b2.Pt();
            m_iter_data->b2_eta = b2.Eta();
            m_iter_data->b2_phi = b2.Phi();
            m_iter_data->b2_mass = b2.M();

            m_iter_data->Hbb_pt = Hbb.Pt();
            m_iter_data->Hbb_eta = Hbb.Eta();
            m_iter_data->Hbb_phi = Hbb.Phi();
            m_iter_data->Hbb_mass = Hbb.M();

            m_iter_data->bjet_resc_fact_1 = c1;
            m_iter_data->bjet_resc_fact_2 = c2;
            m_iter_data->mh = mh;
            m_iter_data->mw = mw;
            m_iter_data->smear_dpx = smear_dpx;
            m_iter_data->smear_dpy = smear_dpy;
            m_iter_data->bjet_resc_dpx = bjet_resc_dpx;
            m_iter_data->bjet_resc_dpy = bjet_resc_dpy;
            m_iter_data->weight = weight;
            m_iter_data->num_sol = num_sol;

            m_recorder.FillTree();
            ResetObject(*m_iter_data);
        }
    }

    Float_t integral = m_res_mass->Integral();
    if (m_res_mass->GetEntries() && integral > 0.0)
    {
        if (m_recorder.ShouldRecord())
        {
            m_recorder.WriteTree(tree_name);
        }

        int binmax = m_res_mass->GetMaximumBin(); 
        res[static_cast<size_t>(EstimOut::mass)] = m_res_mass->GetXaxis()->GetBinCenter(binmax);
        res[static_cast<size_t>(EstimOut::peak_value)] = m_res_mass->GetBinContent(binmax);
        res[static_cast<size_t>(EstimOut::width)] = ComputeWidth(m_res_mass, Q16, Q84);
        res[static_cast<size_t>(EstimOut::integral)] = integral;
        return res;
    }
    return res;
}

OptArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id)
{
    VecLVF_t particles(static_cast<size_t>(ObjDL::count));
    particles[static_cast<size_t>(ObjDL::lep1)] = leptons[static_cast<size_t>(Lep::lep1)];
    particles[static_cast<size_t>(ObjDL::lep2)] = leptons[static_cast<size_t>(Lep::lep2)];
    particles[static_cast<size_t>(ObjDL::met)] = met;

    std::vector<ArrF_t<ESTIM_OUT_SZ>> estimations;
    size_t num_bjets = jets.size() < NUM_BEST_BTAG ? jets.size() : NUM_BEST_BTAG;

    if (m_aggr_mode == AggregationMode::Event)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu", evt_id));
    }
    
    for (size_t bj1_idx = 0; bj1_idx < num_bjets; ++bj1_idx)
    {
        for (size_t bj2_idx = bj1_idx + 1; bj2_idx < num_bjets; ++bj2_idx)
        {
            // order jets such that first b jet has bigger pt and save their p4
            if (jets[bj1_idx].Pt() > jets[bj2_idx].Pt())
            {
                particles[static_cast<size_t>(ObjDL::bj1)] = jets[bj1_idx];
                particles[static_cast<size_t>(ObjDL::bj2)] = jets[bj2_idx];
            }
            else 
            {
                particles[static_cast<size_t>(ObjDL::bj1)] = jets[bj2_idx];
                particles[static_cast<size_t>(ObjDL::bj2)] = jets[bj1_idx];
            }

            TString comb_label = Form("b%zub%zu", bj1_idx, bj2_idx);
            ArrF_t<ESTIM_OUT_SZ> comb_result = EstimateCombination(particles, evt_id, comb_label);

            // success: mass > 0
            if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
            {
                estimations.push_back(comb_result);
            }

            // clear the histogram to be reused 
            if (m_aggr_mode == AggregationMode::Combination)
            {
                ResetHist(m_res_mass);
            }
        }
    }

    // success: at least one combination produced an estimate of X->HH mass
    if (!estimations.empty())
    {
        if (m_aggr_mode == AggregationMode::Combination)
        {
            ResetHist(m_res_mass);
            return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(estimations[0]);
        }
        else if (m_aggr_mode == AggregationMode::Event)
        {
            ArrF_t<ESTIM_OUT_SZ> res{};
            int binmax = m_res_mass->GetMaximumBin(); 
            res[static_cast<size_t>(EstimOut::mass)] = m_res_mass->GetXaxis()->GetBinCenter(binmax);
            res[static_cast<size_t>(EstimOut::peak_value)] = m_res_mass->GetBinContent(binmax);
            res[static_cast<size_t>(EstimOut::width)] = ComputeWidth(m_res_mass, Q16, Q84);
            res[static_cast<size_t>(EstimOut::integral)] = m_res_mass->Integral();
            ResetHist(m_res_mass);
            return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(res);
        }
        else 
        {
            throw std::runtime_error("Unknown strategy to aggregate combination data to estimate event mass");
        }
    }
    ResetHist(m_res_mass);
    return std::nullopt;
}

std::unique_ptr<TTree> EstimatorDoubleLep::MakeTree(TString const& tree_name)
{
    auto tree = std::make_unique<TTree>(tree_name, "DL channel debug tree");
    tree->SetDirectory(nullptr);

    tree->Branch("offshellW_pt", m_iter_data->offshellW_pt, "offshellW_pt[4]/F");
    tree->Branch("offshellW_eta", m_iter_data->offshellW_eta, "offshellW_eta[4]/F");
    tree->Branch("offshellW_phi", m_iter_data->offshellW_phi, "offshellW_phi[4]/F");
    tree->Branch("offshellW_mass", m_iter_data->offshellW_mass, "offshellW_mass[4]/F");

    tree->Branch("onshellW_pt", m_iter_data->onshellW_pt, "onshellW_pt[4]/F");
    tree->Branch("onshellW_eta", m_iter_data->onshellW_eta, "onshellW_eta[4]/F");
    tree->Branch("onshellW_phi", m_iter_data->onshellW_phi, "onshellW_phi[4]/F");
    tree->Branch("onshellW_mass", m_iter_data->onshellW_mass, "onshellW_mass[4]/F");

    tree->Branch("met_corr_pt", &m_iter_data->met_corr_pt, "met_corr_pt/F");
    tree->Branch("met_corr_phi", &m_iter_data->met_corr_phi, "met_corr_phi/F");

    tree->Branch("nu_offshell_pt", m_iter_data->nu_offshell_pt, "nu_offshell_pt[4]/F");
    tree->Branch("nu_offshell_eta", m_iter_data->nu_offshell_eta, "nu_offshell_eta[4]/F");
    tree->Branch("nu_offshell_phi", m_iter_data->nu_offshell_phi, "nu_offshell_phi[4]/F");

    tree->Branch("nu_onshell_pt", m_iter_data->nu_onshell_pt, "nu_onshell_pt[4]/F");
    tree->Branch("nu_onshell_eta", m_iter_data->nu_onshell_eta, "nu_onshell_eta[4]/F");
    tree->Branch("nu_onshell_phi", m_iter_data->nu_onshell_phi, "nu_onshell_phi[4]/F");

    tree->Branch("l_offshell_pt", m_iter_data->l_offshell_pt, "l_offshell_pt[4]/F");
    tree->Branch("l_offshell_eta", m_iter_data->l_offshell_eta, "l_offshell_eta[4]/F");
    tree->Branch("l_offshell_phi", m_iter_data->l_offshell_phi, "l_offshell_phi[4]/F");
    tree->Branch("l_offshell_mass", m_iter_data->l_offshell_mass, "l_offshell_mass[4]/F");

    tree->Branch("l_onshell_pt", m_iter_data->l_onshell_pt, "l_onshell_pt[4]/F");
    tree->Branch("l_onshell_eta", m_iter_data->l_onshell_eta, "l_onshell_eta[4]/F");
    tree->Branch("l_onshell_phi", m_iter_data->l_onshell_phi, "l_onshell_phi[4]/F");
    tree->Branch("l_onshell_mass", m_iter_data->l_onshell_mass, "l_onshell_mass[4]/F");

    tree->Branch("Hww_pt", m_iter_data->Hww_pt, "Hww_pt[4]/F");
    tree->Branch("Hww_eta", m_iter_data->Hww_eta, "Hww_eta[4]/F");
    tree->Branch("Hww_phi", m_iter_data->Hww_phi, "Hww_phi[4]/F");
    tree->Branch("Hww_mass", m_iter_data->Hww_mass, "Hww_mass[4]/F");

    tree->Branch("Xhh_pt", m_iter_data->Xhh_pt, "Xhh_pt[4]/F");
    tree->Branch("Xhh_eta", m_iter_data->Xhh_eta, "Xhh_eta[4]/F");
    tree->Branch("Xhh_phi", m_iter_data->Xhh_phi, "Xhh_phi[4]/F");
    tree->Branch("Xhh_mass", m_iter_data->Xhh_mass, "Xhh_mass[4]/F");

    tree->Branch("bjet_1_pt", &m_iter_data->b1_pt, "bjet_1_pt/F");
    tree->Branch("bjet_1_eta", &m_iter_data->b1_eta, "bjet_1_eta/F");
    tree->Branch("bjet_1_phi", &m_iter_data->b1_phi, "bjet_1_phi/F");
    tree->Branch("bjet_1_mass", &m_iter_data->b1_mass, "bjet_1_mass/F");

    tree->Branch("bjet_2_pt", &m_iter_data->b2_pt, "bjet_2_pt/F");
    tree->Branch("bjet_2_eta", &m_iter_data->b2_eta, "bjet_2_eta/F");
    tree->Branch("bjet_2_phi", &m_iter_data->b2_phi, "bjet_2_phi/F");
    tree->Branch("bjet_2_mass", &m_iter_data->b2_mass, "bjet_2_mass/F");

    tree->Branch("Hbb_pt", &m_iter_data->Hbb_pt, "Hbb_pt/F");
    tree->Branch("Hbb_eta", &m_iter_data->Hbb_eta, "Hbb_eta/F");
    tree->Branch("Hbb_phi", &m_iter_data->Hbb_phi, "Hbb_phi/F");
    tree->Branch("Hbb_mass", &m_iter_data->Hbb_mass, "Hbb_mass/F");
    
    tree->Branch("bjet_resc_fact_1", &m_iter_data->bjet_resc_fact_1, "bjet_resc_fact_1/F");
    tree->Branch("bjet_resc_fact_2", &m_iter_data->bjet_resc_fact_2, "bjet_resc_fact_2/F");
    tree->Branch("mh", &m_iter_data->mh, "mh/F");
    tree->Branch("mw", &m_iter_data->mw, "mw/F");
    tree->Branch("smear_dpx", &m_iter_data->smear_dpx, "smear_dpx/F");
    tree->Branch("smear_dpy", &m_iter_data->smear_dpy, "smear_dpy/F");
    tree->Branch("bjet_resc_dpx", &m_iter_data->bjet_resc_dpx, "bjet_resc_dpx/F");
    tree->Branch("bjet_resc_dpy", &m_iter_data->bjet_resc_dpy, "bjet_resc_dpy/F");
    tree->Branch("mass", m_iter_data->mass, "mass[4]/F");
    tree->Branch("weight", &m_iter_data->weight, "weight/F");
    tree->Branch("num_sol", &m_iter_data->num_sol, "num_sol/I");
    return tree;
}