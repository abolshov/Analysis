#include "EstimatorSingleLep.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_set>

#include "TVector2.h"
#include "Math/GenVector/VectorUtil.h" // DeltaPhi
using ROOT::Math::VectorUtil::DeltaPhi;

#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"

EstimatorSingleLep::EstimatorSingleLep(TString const& pdf_file_name, AggregationMode aggr_mode)
:   EstimatorBase(aggr_mode)
,   m_iter_data(std::make_unique<IterData>())
{
    m_pdf_1d.resize(pdf1d_sl_names.size());
    m_pdf_2d.resize(pdf2d_sl_names.size());

    TFile* pf = TFile::Open(pdf_file_name);
    Get1dPDFs(pf, m_pdf_1d, Channel::SL);
    Get2dPDFs(pf, m_pdf_2d, Channel::SL);
    pf->Close();
}

ArrF_t<ESTIM_OUT_SZ> EstimatorSingleLep::EstimateCombSlim(VecLVF_t const& particles, std::vector<Float_t> const& jet_res, ULong64_t evt_id, JetComb const& comb)
{
    ArrF_t<ESTIM_OUT_SZ> res{};
    std::fill(res.begin(), res.end(), -1.0f);

    LorentzVectorF_t const& bj1 = particles[static_cast<size_t>(ObjSL::bj1)];
    LorentzVectorF_t const& bj2 = particles[static_cast<size_t>(ObjSL::bj2)];
    LorentzVectorF_t const& lj1 = particles[static_cast<size_t>(ObjSL::lj1)];
    LorentzVectorF_t const& lj2 = particles[static_cast<size_t>(ObjSL::lj2)];
    LorentzVectorF_t const& lep = particles[static_cast<size_t>(ObjSL::lep)];
    LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjSL::met)];

    UHist_t<TH1F>& pdf_b1 = m_pdf_1d[static_cast<size_t>(PDF1_sl::b1)];
    UHist_t<TH1F>& pdf_q1 = m_pdf_1d[static_cast<size_t>(PDF1_sl::q1)];
    UHist_t<TH2F>& pdf_mw1mw2 = m_pdf_2d[static_cast<size_t>(PDF2_sl::mw1mw2)];

    TString label = comb.ToString(Channel::SL);
    Float_t proba = comb.GetProbability();
    Float_t mh = m_prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);

    if (m_aggr_mode == AggregationMode::Combination)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt_id, label.Data()));
    }
    
    [[maybe_unused]] TString tree_name = Form("evt_%llu_%s", evt_id, label.Data());
    if (m_recorder.ShouldRecord())
    {
        m_recorder.ResetTree(MakeTree(tree_name));
    }

    [[maybe_unused]] int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        Float_t unclust_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t unclust_dpy = m_prg->Gaus(0.0, MET_SIGMA);

        auto bresc = ComputeJetResc(bj1, bj2, pdf_b1, mh);
        if (!bresc.has_value())
        {
            ++failed_iter;
            continue;
        }
        auto [c1, c2] = bresc.value();

        LorentzVectorF_t b1 = bj1;
        LorentzVectorF_t b2 = bj2;
        b1 *= c1;
        b2 *= c2;

        LorentzVectorF_t Hbb = b1 + b2;

        Float_t bjet_resc_dpx = -1.0*(c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px();
        Float_t bjet_resc_dpy = -1.0*(c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py();

        Double_t mw1 = 1.0;
        Double_t mw2 = 1.0;
        pdf_mw1mw2->GetRandom2(mw1, mw2, m_prg.get());

        std::vector<Float_t> masses;
        std::array<LorentzVectorF_t, CONTROL> j1{};
        std::array<LorentzVectorF_t, CONTROL> j2{};
        std::array<LorentzVectorF_t, CONTROL> met_corr{};
        std::array<LorentzVectorF_t, CONTROL> lepW{};
        std::array<LorentzVectorF_t, CONTROL> hadW{};
        std::array<LorentzVectorF_t, CONTROL> Hww{};
        std::array<LorentzVectorF_t, CONTROL> Xhh{};
        std::array<LorentzVectorF_t, CONTROL> nu{};
        std::array<Float_t, CONTROL> c3{};
        std::array<Float_t, CONTROL> c4{};
        std::array<Float_t, CONTROL> ljet_resc_dpx{};
        std::array<Float_t, CONTROL> ljet_resc_dpy{};
        std::array<Float_t, CONTROL> mass{};
        std::array<Bool_t, CONTROL> correct_hww_mass{};
        for (int control = 0; control < CONTROL; ++control)
        {
            bool lepW_onshell = control / 2;
            bool add_deta = control % 2;

            Float_t mWlep = lepW_onshell ? mw1 : mw2;
            Float_t mWhad = lepW_onshell ? mw2 : mw1;

            j1[control] = lj1;
            j2[control] = lj2;
            auto lresc = ComputeJetResc(j1[control], j2[control], pdf_q1, mWhad);
            if (!lresc.has_value())
            {
                ++failed_iter;
                continue;
            }

            std::tie(c3[control], c4[control]) = lresc.value();
            j1[control] *= c3[control];
            j2[control] *= c4[control];
            
            ljet_resc_dpx[control] = -1.0*(c3[control] - 1)*lj1.Px() - (c4[control] - 1)*lj2.Px();
            ljet_resc_dpy[control] = -1.0*(c3[control] - 1)*lj1.Py() - (c4[control] - 1)*lj2.Py();

            Float_t met_corr_px = met.Px() + bjet_resc_dpx + ljet_resc_dpx[control] + unclust_dpx;
            Float_t met_corr_py = met.Py() + bjet_resc_dpy + ljet_resc_dpy[control] + unclust_dpy;

            Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
            Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
            met_corr[control] = LorentzVectorF_t(met_corr_pt, 0.0, met_corr_phi, 0.0);
            
            auto opt_nu = NuFromW(lep, met_corr[control], add_deta, mWlep);
            if (opt_nu)
            {
                nu[control] = opt_nu.value();
                lepW[control] = nu[control] + lep;
                hadW[control] = j1[control] + j2[control];
                Hww[control] = lepW[control] + hadW[control];
                Xhh[control] = Hww[control] + Hbb;
                mass[control] = Xhh[control].M();
                correct_hww_mass[control] = (std::abs(mh - Hww[control].M()) < 1.0);
                if (!correct_hww_mass[control])
                {
                    continue;
                }
                masses.push_back(mass[control]);
            } 
            else 
            {
                continue;
            }
        }

        if (masses.empty())
        {
            ++failed_iter;
            continue;
        }

        Int_t num_sol = masses.size();
        Float_t weight = 1.0/num_sol;
        for (auto mass: masses)
        {
            m_res_mass->Fill(mass, proba*weight);
        }
        
        if (m_recorder.ShouldRecord())
        {
            // here iteration was succusseful
            // assign local vaiables to members of IterData and dump them in a tree
            // after writing, reset all members of IterData to avoid writing junk values from previous iterations 
            for (int i = 0; i < CONTROL; ++i)
            {
                m_iter_data->j1_pt[i] = j1[i].Pt();
                m_iter_data->j1_eta[i] = j1[i].Eta();
                m_iter_data->j1_phi[i] = j1[i].Phi();
                m_iter_data->j1_mass[i] = j1[i].M();

                m_iter_data->j2_pt[i] = j2[i].Pt();
                m_iter_data->j2_eta[i] = j2[i].Eta();
                m_iter_data->j2_phi[i] = j2[i].Phi();
                m_iter_data->j2_mass[i] = j2[i].M();

                m_iter_data->lepW_pt[i] = lepW[i].Pt();
                m_iter_data->lepW_eta[i] = lepW[i].Eta();
                m_iter_data->lepW_phi[i] = lepW[i].Phi();
                m_iter_data->lepW_mass[i] = lepW[i].M();

                m_iter_data->hadW_pt[i] = hadW[i].Pt();
                m_iter_data->hadW_eta[i] = hadW[i].Eta();
                m_iter_data->hadW_phi[i] = hadW[i].Phi();
                m_iter_data->hadW_mass[i] = hadW[i].M();

                m_iter_data->Hww_pt[i] = Hww[i].Pt();
                m_iter_data->Hww_eta[i] = Hww[i].Eta();
                m_iter_data->Hww_phi[i] = Hww[i].Phi();
                m_iter_data->Hww_mass[i] = Hww[i].M();

                m_iter_data->nu_pt[i] = nu[i].Pt();
                m_iter_data->nu_eta[i] = nu[i].Eta();
                m_iter_data->nu_phi[i] = nu[i].Phi();

                m_iter_data->Xhh_pt[i] = Xhh[i].Pt();
                m_iter_data->Xhh_eta[i] = Xhh[i].Eta();
                m_iter_data->Xhh_phi[i] = Xhh[i].Phi();
                m_iter_data->Xhh_mass[i] = Xhh[i].M();

                m_iter_data->met_corr_pt[i] = met_corr[i].Pt();
                m_iter_data->met_corr_phi[i] = met_corr[i].Phi();

                m_iter_data->mass[i] = mass[i];

                m_iter_data->ljet_resc_fact_1[i] = c3[i];
                m_iter_data->ljet_resc_fact_2[i] = c4[i];

                m_iter_data->ljet_resc_dpx[i] = ljet_resc_dpx[i];
                m_iter_data->ljet_resc_dpy[i] = ljet_resc_dpy[i];

                m_iter_data->correct_hww_mass[i] = correct_hww_mass[i];
            }

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
            m_iter_data->mw1 = mw1;
            m_iter_data->mw2 = mw2;
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

    // combination data is returned in any case for further analysis
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

OptArrF_t<ESTIM_OUT_SZ> EstimatorSingleLep::EstimateMass(Event const& event)
{
    using JetIdxPair_t = std::pair<size_t, size_t>; 

    std::vector<LorentzVectorF_t> particles(static_cast<size_t>(ObjSL::count));
    particles[static_cast<size_t>(ObjSL::lep)] = LorentzVectorF_t(event.reco_lep_pt[static_cast<size_t>(Lep::lep1)], 
                                                                  event.reco_lep_eta[static_cast<size_t>(Lep::lep1)],
                                                                  event.reco_lep_phi[static_cast<size_t>(Lep::lep1)], 
                                                                  event.reco_lep_mass[static_cast<size_t>(Lep::lep1)]);
    particles[static_cast<size_t>(ObjSL::met)] = LorentzVectorF_t(event.reco_met_pt, 0.0f, event.reco_met_phi, 0.0f);
    
    if (m_aggr_mode == AggregationMode::Event)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu", event.event_id));
    }

    size_t n_jets = static_cast<size_t>(event.n_reco_jet);
    size_t num_bjets = n_jets < NUM_BEST_BTAG ? n_jets : NUM_BEST_BTAG;
    std::map<JetIdxPair_t, std::vector<JetIdxPair_t>> jet_comb_candidates;
    for (size_t i = 0; i < num_bjets; ++i)
    {
        for (size_t j = i + 1; j < num_bjets; ++j)
        {
            jet_comb_candidates.insert({{i, j}, std::vector<JetIdxPair_t>()});
        }   
    }

    // fill possible light jet indices
    // 2 jets are used to form a bjet pair => n - 2 jets can be light jet candidates
    size_t num_light_jets = USE_QVG_SELECTION ? std::min(n_jets - 2, NUM_BEST_QVG) : n_jets - 2;
    for (auto& [bjet_pair, light_jet_idx_pairs]: jet_comb_candidates)
    {
        auto [bj1_idx, bj2_idx] = bjet_pair;
        std::vector<size_t> light_jet_indices;
        for (size_t lj_idx = 0; lj_idx < n_jets; ++lj_idx)
        {
            if (lj_idx == bj1_idx || lj_idx == bj2_idx)
            {
                continue;
            }
            light_jet_indices.push_back(lj_idx);
        }
        
        // sort light jet indices by btagPNetQvG to easily select n best
        if (USE_QVG_SELECTION)
        {
            std::array<Float_t, MAX_RECO_JET> const& qvg_scores = event.reco_jet_qvg;
            auto ScoreCmp = [&qvg_scores](size_t i, size_t j)
            {
                return qvg_scores[i] > qvg_scores[j];
            };
            std::sort(light_jet_indices.begin(), light_jet_indices.end(), ScoreCmp);
        }

        // form all possible pairs of light jets
        for (size_t i = 0; i < num_light_jets; ++i)
        {
            for (size_t j = i + 1; j < num_light_jets; ++j)
            {
                light_jet_idx_pairs.emplace_back(light_jet_indices[i], light_jet_indices[j]);
            }   
        }
    }

    // loop over combination and estimate mass for each
    std::vector<ArrF_t<ESTIM_OUT_SZ>> estimations;
    for (auto const& [bjet_idx_pair, light_jet_idx_pairs]: jet_comb_candidates)
    {
        auto [bj1_idx, bj2_idx] = bjet_idx_pair;
        if (event.reco_jet_pt[bj1_idx] > event.reco_jet_pt[bj2_idx])
        {
            particles[static_cast<size_t>(ObjSL::bj1)] = LorentzVectorF_t(event.reco_jet_pt[bj1_idx], 
                                                                          event.reco_jet_eta[bj1_idx],
                                                                          event.reco_jet_phi[bj1_idx], 
                                                                          event.reco_jet_mass[bj1_idx]);
            particles[static_cast<size_t>(ObjSL::bj2)] = LorentzVectorF_t(event.reco_jet_pt[bj2_idx], 
                                                                          event.reco_jet_eta[bj2_idx],
                                                                          event.reco_jet_phi[bj2_idx], 
                                                                          event.reco_jet_mass[bj2_idx]);
        }
        else 
        {
            particles[static_cast<size_t>(ObjSL::bj1)] = LorentzVectorF_t(event.reco_jet_pt[bj2_idx], 
                                                                          event.reco_jet_eta[bj2_idx],
                                                                          event.reco_jet_phi[bj2_idx], 
                                                                          event.reco_jet_mass[bj2_idx]);
            particles[static_cast<size_t>(ObjSL::bj2)] = LorentzVectorF_t(event.reco_jet_pt[bj1_idx], 
                                                                          event.reco_jet_eta[bj1_idx],
                                                                          event.reco_jet_phi[bj1_idx], 
                                                                          event.reco_jet_mass[bj1_idx]);
        }

        for (auto const& [lj1_idx, lj2_idx]: light_jet_idx_pairs)
        {
            std::vector<Float_t> other_jet_resolutions{}; // dummy for now

            if (event.reco_jet_pt[lj1_idx] > event.reco_jet_pt[lj2_idx])
            {
                particles[static_cast<size_t>(ObjSL::lj1)] = LorentzVectorF_t(event.reco_jet_pt[lj1_idx], 
                                                                              event.reco_jet_eta[lj1_idx],
                                                                              event.reco_jet_phi[lj1_idx], 
                                                                              event.reco_jet_mass[lj1_idx]);
                particles[static_cast<size_t>(ObjSL::lj2)] = LorentzVectorF_t(event.reco_jet_pt[lj2_idx], 
                                                                              event.reco_jet_eta[lj2_idx],
                                                                              event.reco_jet_phi[lj2_idx], 
                                                                              event.reco_jet_mass[lj2_idx]);
            }
            else 
            {
                particles[static_cast<size_t>(ObjSL::lj1)] = LorentzVectorF_t(event.reco_jet_pt[lj2_idx], 
                                                                              event.reco_jet_eta[lj2_idx],
                                                                              event.reco_jet_phi[lj2_idx], 
                                                                              event.reco_jet_mass[lj2_idx]);
                particles[static_cast<size_t>(ObjSL::lj2)] = LorentzVectorF_t(event.reco_jet_pt[lj1_idx], 
                                                                              event.reco_jet_eta[lj1_idx],
                                                                              event.reco_jet_phi[lj1_idx], 
                                                                              event.reco_jet_mass[lj1_idx]);
            }

            JetComb comb(bj1_idx, bj2_idx, lj1_idx, lj2_idx, event.reco_jet_btag, event.reco_jet_qvg);
            ArrF_t<ESTIM_OUT_SZ> comb_result = EstimateCombSlim(particles, other_jet_resolutions, event.event_id, comb);

            // success: mass > 0
            if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
            {
                estimations.push_back(comb_result);
            }
        }
    }

    if (!estimations.empty())
    {
        if (m_aggr_mode == AggregationMode::Combination)
        {
            // pick combination with largest peak_value
            auto PeakComparator = [](ArrF_t<ESTIM_OUT_SZ> const& c1, ArrF_t<ESTIM_OUT_SZ> const& c2)
            {
                return c1[static_cast<size_t>(EstimOut::peak_value)] < c2[static_cast<size_t>(EstimOut::peak_value)];
            };
            auto it = std::max_element(estimations.begin(), estimations.end(), PeakComparator);
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

std::unique_ptr<TTree> EstimatorSingleLep::MakeTree(TString const& tree_name)
{
    auto tree = std::make_unique<TTree>(tree_name, "SL channel debug tree");
    tree->SetDirectory(nullptr);

    tree->Branch("light_jet_1_pt", m_iter_data->j1_pt, "light_jet_1_pt[4]/F");
    tree->Branch("light_jet_1_eta", m_iter_data->j1_eta, "light_jet_1_eta[4]/F");
    tree->Branch("light_jet_1_phi", m_iter_data->j1_phi, "light_jet_1_phi[4]/F");
    tree->Branch("light_jet_1_mass", m_iter_data->j1_mass, "light_jet_1_mass[4]/F");

    tree->Branch("light_jet_2_pt", m_iter_data->j2_pt, "light_jet_2_pt[4]/F");
    tree->Branch("light_jet_2_eta", m_iter_data->j2_eta, "light_jet_2_eta[4]/F");
    tree->Branch("light_jet_2_phi", m_iter_data->j2_phi, "light_jet_2_phi[4]/F");
    tree->Branch("light_jet_2_mass", m_iter_data->j2_mass, "light_jet_2_mass[4]/F");

    tree->Branch("met_corr_pt", m_iter_data->met_corr_pt, "met_corr_pt[4]/F");
    tree->Branch("met_corr_phi", m_iter_data->met_corr_phi, "met_corr_phi[4]/F");

    tree->Branch("nu_pt", m_iter_data->nu_pt, "nu_pt[4]/F");
    tree->Branch("nu_eta", m_iter_data->nu_eta, "nu_eta[4]/F");
    tree->Branch("nu_phi", m_iter_data->nu_phi, "nu_phi[4]/F");

    tree->Branch("lepW_pt", m_iter_data->lepW_pt, "lepW_pt[4]/F");
    tree->Branch("lepW_eta", m_iter_data->lepW_eta, "lepW_eta[4]/F");
    tree->Branch("lepW_phi", m_iter_data->lepW_phi, "lepW_phi[4]/F");
    tree->Branch("lepW_mass", m_iter_data->lepW_mass, "lepW_mass[4]/F");

    tree->Branch("hadW_pt", m_iter_data->hadW_pt, "hadW_pt[4]/F");
    tree->Branch("hadW_eta", m_iter_data->hadW_eta, "hadW_eta[4]/F");
    tree->Branch("hadW_phi", m_iter_data->hadW_phi, "hadW_phi[4]/F");
    tree->Branch("hadW_mass", m_iter_data->hadW_mass, "hadW_mass[4]/F");

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
    tree->Branch("mw1", &m_iter_data->mw1, "mw1/F");
    tree->Branch("mw2", &m_iter_data->mw2, "mw2/F");
    tree->Branch("unclust_dpx", &m_iter_data->unclust_dpx, "unclust_dpx/F");
    tree->Branch("unclust_dpy", &m_iter_data->unclust_dpy, "unclust_dpy/F");
    tree->Branch("bjet_resc_dpx", &m_iter_data->bjet_resc_dpx, "bjet_resc_dpx/F");
    tree->Branch("bjet_resc_dpy", &m_iter_data->bjet_resc_dpy, "bjet_resc_dpy/F");
    tree->Branch("ljet_resc_fact_1", m_iter_data->ljet_resc_fact_1, "ljet_resc_fact_1[4]/F");
    tree->Branch("ljet_resc_fact_2", m_iter_data->ljet_resc_fact_2, "ljet_resc_fact_2[4]/F");
    tree->Branch("ljet_resc_dpx", m_iter_data->ljet_resc_dpx, "ljet_resc_dpx[4]/F");
    tree->Branch("ljet_resc_dpy", m_iter_data->ljet_resc_dpy, "ljet_resc_dpy[4]/F");
    tree->Branch("mass", m_iter_data->mass, "mass[4]/F");
    tree->Branch("weight", &m_iter_data->weight, "weight/F");
    tree->Branch("num_sol", &m_iter_data->num_sol, "num_sol/I");
    tree->Branch("correct_hww_mass", m_iter_data->correct_hww_mass, "correct_hww_mass[4]/B");
    return tree;
}