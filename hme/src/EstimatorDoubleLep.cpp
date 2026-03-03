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

ArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateCombSlim(VecLVF_t const& particles, 
                                                          [[maybe_unused]] std::vector<Float_t> const& jet_res, 
                                                          ULong64_t evt_id, 
                                                          JetComb const& comb)
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

    TString label = comb.ToString(Channel::DL);
    Float_t proba = comb.GetProbability();

    if (m_aggr_mode == AggregationMode::Combination)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt_id, label.Data()));
    }

    [[maybe_unused]] TString tree_name = Form("evt_%llu_%s", evt_id, label.Data());
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
        Float_t unclust_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t unclust_dpy = m_prg->Gaus(0.0, MET_SIGMA);

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

        Float_t met_corr_px = met.Px() + bjet_resc_dpx + unclust_dpx;
        Float_t met_corr_py = met.Py() + bjet_resc_dpy + unclust_dpy;

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
            m_res_mass->Fill(est, proba*weight);
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
            m_iter_data->unclust_dpx = unclust_dpx;
            m_iter_data->unclust_dpy = unclust_dpy;
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

        if (m_aggr_mode == AggregationMode::Combination)
        {
            ResetHist(m_res_mass);
        }

        return res;
    }

    if (m_aggr_mode == AggregationMode::Combination)
    {
        ResetHist(m_res_mass);
    }
    return res;
}

// only implements logic of using slim (ak4) jets
OptArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateMass(Event const& event)
{
    std::vector<LorentzVectorF_t> particles(static_cast<size_t>(ObjDL::count));
    particles[static_cast<size_t>(ObjDL::lep1)] = LorentzVectorF_t(event.reco_lep_pt[static_cast<size_t>(Lep::lep1)], 
                                                                   event.reco_lep_eta[static_cast<size_t>(Lep::lep1)],
                                                                   event.reco_lep_phi[static_cast<size_t>(Lep::lep1)], 
                                                                   event.reco_lep_mass[static_cast<size_t>(Lep::lep1)]);
    particles[static_cast<size_t>(ObjDL::lep2)] = LorentzVectorF_t(event.reco_lep_pt[static_cast<size_t>(Lep::lep2)], 
                                                                   event.reco_lep_eta[static_cast<size_t>(Lep::lep2)],
                                                                   event.reco_lep_phi[static_cast<size_t>(Lep::lep2)], 
                                                                   event.reco_lep_mass[static_cast<size_t>(Lep::lep2)]);
    particles[static_cast<size_t>(ObjDL::met)] = LorentzVectorF_t(event.reco_met_pt, 0.0f, event.reco_met_phi, 0.0f);
    
    if (m_aggr_mode == AggregationMode::Event)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu", event.event_id));
    }

    // selection of b jets
    std::vector<std::pair<size_t, size_t>> bjet_pair_indices;
    size_t num_bjets = static_cast<size_t>(event.n_reco_jet) < NUM_BEST_BTAG ? static_cast<size_t>(event.n_reco_jet) : NUM_BEST_BTAG;
    for (size_t i = 0; i < num_bjets; ++i)
    {
        for (size_t j = i + 1; j < num_bjets; ++j)
        {
            bjet_pair_indices.emplace_back(i, j);
        }   
    }

    // now loop over all saved jet pairs
    // construct their p4 and push back it to particles
    // and estimate each combination
    std::vector<ArrF_t<ESTIM_OUT_SZ>> estimations;
    for (auto const& [bj1_idx, bj2_idx]: bjet_pair_indices)
    {
        // PNet resolutions of all jets but currently selected b jets
        // will be needed in fututre when I add contributions to met corrections from other jets
        // right now just a placeholder
        std::vector<Float_t> other_jet_resolutions{};

        // order jets such that first b jet has bigger pt and save their p4
        if (event.reco_jet_pt[bj1_idx] > event.reco_jet_pt[bj2_idx])
        {
            particles[static_cast<size_t>(ObjDL::bj1)] = LorentzVectorF_t(event.reco_jet_pt[bj1_idx], 
                                                                          event.reco_jet_eta[bj1_idx],
                                                                          event.reco_jet_phi[bj1_idx], 
                                                                          event.reco_jet_mass[bj1_idx]);
            particles[static_cast<size_t>(ObjDL::bj2)] = LorentzVectorF_t(event.reco_jet_pt[bj2_idx], 
                                                                          event.reco_jet_eta[bj2_idx],
                                                                          event.reco_jet_phi[bj2_idx], 
                                                                          event.reco_jet_mass[bj2_idx]);
        }
        else 
        {
            particles[static_cast<size_t>(ObjDL::bj1)] = LorentzVectorF_t(event.reco_jet_pt[bj2_idx], 
                                                                          event.reco_jet_eta[bj2_idx],
                                                                          event.reco_jet_phi[bj2_idx], 
                                                                          event.reco_jet_mass[bj2_idx]);
            particles[static_cast<size_t>(ObjDL::bj2)] = LorentzVectorF_t(event.reco_jet_pt[bj1_idx], 
                                                                          event.reco_jet_eta[bj1_idx],
                                                                          event.reco_jet_phi[bj1_idx], 
                                                                          event.reco_jet_mass[bj1_idx]);
        }
        
        JetComb comb(bj1_idx, bj2_idx, event.reco_jet_btag);
        ArrF_t<ESTIM_OUT_SZ> comb_result = EstimateCombSlim(particles, other_jet_resolutions, event.event_id, comb);

        // success: mass > 0
        if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
        {
            estimations.push_back(comb_result);
        }
    }

    if (!estimations.empty())
    {
        if (m_aggr_mode == AggregationMode::Combination)
        {
            // pick combination with largest integral
            auto IntegralComparator = [](ArrF_t<ESTIM_OUT_SZ> const& c1, ArrF_t<ESTIM_OUT_SZ> const& c2)
            {
                return c1[static_cast<size_t>(EstimOut::integral)] < c2[static_cast<size_t>(EstimOut::integral)];
            };
            auto it = std::max_element(estimations.begin(), estimations.end(), IntegralComparator);
            size_t best_est_idx = it - estimations.begin();
            return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(estimations[best_est_idx]);
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

// OptArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateMass(Event const& event)
// {
//     std::vector<LorentzVectorF_t> particles(static_cast<size_t>(ObjDL::count));
//     particles[static_cast<size_t>(ObjDL::lep1)] = LorentzVectorF_t(event.reco_lep_pt[static_cast<size_t>(Lep::lep1)], 
//                                                                    event.reco_lep_eta[static_cast<size_t>(Lep::lep1)],
//                                                                    event.reco_lep_phi[static_cast<size_t>(Lep::lep1)], 
//                                                                    event.reco_lep_mass[static_cast<size_t>(Lep::lep1)]);
//     particles[static_cast<size_t>(ObjDL::lep2)] = LorentzVectorF_t(event.reco_lep_pt[static_cast<size_t>(Lep::lep2)], 
//                                                                    event.reco_lep_eta[static_cast<size_t>(Lep::lep2)],
//                                                                    event.reco_lep_phi[static_cast<size_t>(Lep::lep2)], 
//                                                                    event.reco_lep_mass[static_cast<size_t>(Lep::lep2)]);
//     particles[static_cast<size_t>(ObjDL::met)] = LorentzVectorF_t(event.reco_met_pt, 0.0f, event.reco_met_phi, 0.0f);
    
//     if (m_aggr_mode == AggregationMode::Event)
//     {
//         m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu", event.event_id));
//     }

//     size_t num_fatjets = event.n_reco_fatjet;
//     if (num_fatjets == 0)
//     {
//         return std::nullopt;
//     }

//     // select fatjets to consider here
//     // for now only choose the one with highest btag
//     std::vector<ArrF_t<ESTIM_OUT_SZ>> estimations;
//     auto it = std::max_element(event.reco_fatjet_btag.begin(), event.reco_fatjet_btag.end());
//     size_t fatjet_idx = it - event.reco_fatjet_btag.begin();
//     particles[static_cast<size_t>(ObjDL::fatbb)] = LorentzVectorF_t(event.reco_fatjet_pt[fatjet_idx], 
//                                                                     event.reco_fatjet_eta[fatjet_idx],
//                                                                     event.reco_fatjet_phi[fatjet_idx], 
//                                                                     event.reco_fatjet_mass[fatjet_idx]);

//     JetComb comb{};
//     std::vector<Float_t> other_jet_resolutions{};
//     ArrF_t<ESTIM_OUT_SZ> comb_result = EstimateCombFat(particles, other_jet_resolutions, event.event_id, comb);
//     if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
//     {
//         estimations.push_back(comb_result);
//     }

//     if (!estimations.empty())
//     {
//         if (m_aggr_mode == AggregationMode::Combination)
//         {
//             // pick combination with largest integral
//             auto IntegralComparator = [](ArrF_t<ESTIM_OUT_SZ> const& c1, ArrF_t<ESTIM_OUT_SZ> const& c2)
//             {
//                 return c1[static_cast<size_t>(EstimOut::integral)] < c2[static_cast<size_t>(EstimOut::integral)];
//             };
//             auto it = std::max_element(estimations.begin(), estimations.end(), IntegralComparator);
//             size_t best_est_idx = it - estimations.begin();
//             return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(estimations[best_est_idx]);
//         }
//         else if (m_aggr_mode == AggregationMode::Event)
//         {
//             ArrF_t<ESTIM_OUT_SZ> res{};
//             int binmax = m_res_mass->GetMaximumBin(); 
//             res[static_cast<size_t>(EstimOut::mass)] = m_res_mass->GetXaxis()->GetBinCenter(binmax);
//             res[static_cast<size_t>(EstimOut::peak_value)] = m_res_mass->GetBinContent(binmax);
//             res[static_cast<size_t>(EstimOut::width)] = ComputeWidth(m_res_mass, Q16, Q84);
//             res[static_cast<size_t>(EstimOut::integral)] = m_res_mass->Integral();
//             ResetHist(m_res_mass);
//             return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(res);
//         }
//         else 
//         {
//             throw std::runtime_error("Unknown strategy to aggregate combination data to estimate event mass");
//         }
//     }
//     ResetHist(m_res_mass);
//     return std::nullopt;
// }

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

    tree->Branch("fatjet_pt", &m_iter_data->fatjet_pt, "fatjet_pt/F");
    tree->Branch("fatjet_eta", &m_iter_data->fatjet_eta, "fatjet_eta/F");
    tree->Branch("fatjet_phi", &m_iter_data->fatjet_phi, "fatjet_phi/F");
    tree->Branch("fatjet_mass", &m_iter_data->fatjet_mass, "fatjet_mass/F");

    tree->Branch("Hbb_pt", &m_iter_data->Hbb_pt, "Hbb_pt/F");
    tree->Branch("Hbb_eta", &m_iter_data->Hbb_eta, "Hbb_eta/F");
    tree->Branch("Hbb_phi", &m_iter_data->Hbb_phi, "Hbb_phi/F");
    tree->Branch("Hbb_mass", &m_iter_data->Hbb_mass, "Hbb_mass/F");
    
    tree->Branch("bjet_resc_fact_1", &m_iter_data->bjet_resc_fact_1, "bjet_resc_fact_1/F");
    tree->Branch("bjet_resc_fact_2", &m_iter_data->bjet_resc_fact_2, "bjet_resc_fact_2/F");
    tree->Branch("fatjet_resc_fact", &m_iter_data->fatjet_resc_fact, "fatjet_resc_fact/F");
    tree->Branch("mh", &m_iter_data->mh, "mh/F");
    tree->Branch("mw", &m_iter_data->mw, "mw/F");
    tree->Branch("unclust_dpx", &m_iter_data->unclust_dpx, "unclust_dpx/F");
    tree->Branch("unclust_dpy", &m_iter_data->unclust_dpy, "unclust_dpy/F");
    tree->Branch("bjet_resc_dpx", &m_iter_data->bjet_resc_dpx, "bjet_resc_dpx/F");
    tree->Branch("bjet_resc_dpy", &m_iter_data->bjet_resc_dpy, "bjet_resc_dpy/F");
    tree->Branch("mass", m_iter_data->mass, "mass[4]/F");
    tree->Branch("weight", &m_iter_data->weight, "weight/F");
    tree->Branch("num_sol", &m_iter_data->num_sol, "num_sol/I");
    return tree;
}

ArrF_t<ESTIM_OUT_SZ> EstimatorDoubleLep::EstimateCombFat(VecLVF_t const& particles, 
                                                         [[maybe_unused]] std::vector<Float_t> const& jet_res,
                                                         ULong64_t evt_id, 
                                                         [[maybe_unused]] JetComb const& comb)
{
    ArrF_t<ESTIM_OUT_SZ> res{};
    std::fill(res.begin(), res.end(), -1.0f);

    LorentzVectorF_t const& fatjet = particles[static_cast<size_t>(ObjDL::fatbb)];
    LorentzVectorF_t const& lep1 = particles[static_cast<size_t>(ObjDL::lep1)];
    LorentzVectorF_t const& lep2 = particles[static_cast<size_t>(ObjDL::lep2)];
    LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjDL::met)];

    UHist_t<TH1F>& pdf_fatbb = m_pdf_1d[static_cast<size_t>(PDF1_dl::fatbb)];
    UHist_t<TH1F>& pdf_mw_onshell = m_pdf_1d[static_cast<size_t>(PDF1_dl::mw_onshell)];

    TString label = comb.ToString(Channel::DL);
    Float_t proba = comb.GetProbability();

    // this has to be fixed: in case of fatjet JetComb is irrelevant
    if (m_aggr_mode == AggregationMode::Combination)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt_id, label.Data()));
    }

    [[maybe_unused]] TString tree_name = Form("evt_%llu_%s", evt_id, label.Data());
    if (m_recorder.ShouldRecord())
    {
        m_recorder.ResetTree(MakeTree(tree_name));
    }

    for (int i = 0; i < N_ITER; ++i)
    {
        // iter loop
        Float_t c_fat = pdf_fatbb->GetRandom(m_prg.get());
        Float_t mh = m_prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);
        // Float_t c_fat = mh/fatjet.M();
        LorentzVectorF_t Hbb = c_fat*fatjet;
        if (std::abs(Hbb.M() - mh) > 1.0)
        {
            continue;
        }

        // generate other random stuff after checking fatbb correction validity to save time in case it fails
        Float_t eta_gen = m_prg->Uniform(-6, 6);
        Float_t phi_gen = m_prg->Uniform(-3.1415926, 3.1415926);
        Float_t mw = pdf_mw_onshell->GetRandom(m_prg.get());
        Float_t unclust_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t unclust_dpy = m_prg->Gaus(0.0, MET_SIGMA);

        Float_t bjet_resc_dpx = -1.0*(c_fat - 1)*fatjet.Px();
        Float_t bjet_resc_dpy = -1.0*(c_fat - 1)*fatjet.Py();

        Float_t met_corr_px = met.Px() + bjet_resc_dpx + unclust_dpx;
        Float_t met_corr_py = met.Py() + bjet_resc_dpy + unclust_dpy;

        Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        LorentzVectorF_t met_corr(met_corr_pt, 0.0, met_corr_phi, 0.0);

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
            m_res_mass->Fill(est, proba*weight);
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

            m_iter_data->fatjet_pt = fatjet.Pt();
            m_iter_data->fatjet_eta = fatjet.Eta();
            m_iter_data->fatjet_phi = fatjet.Phi();
            m_iter_data->fatjet_mass = fatjet.M();

            m_iter_data->Hbb_pt = Hbb.Pt();
            m_iter_data->Hbb_eta = Hbb.Eta();
            m_iter_data->Hbb_phi = Hbb.Phi();
            m_iter_data->Hbb_mass = Hbb.M();

            m_iter_data->fatjet_resc_fact = c_fat;
            m_iter_data->mh = mh;
            m_iter_data->mw = mw;
            m_iter_data->unclust_dpx = unclust_dpx;
            m_iter_data->unclust_dpy = unclust_dpy;
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

        if (m_aggr_mode == AggregationMode::Combination)
        {
            ResetHist(m_res_mass);
        }

        return res;
    }

    if (m_aggr_mode == AggregationMode::Combination)
    {
        ResetHist(m_res_mass);
    }

    return res;
}