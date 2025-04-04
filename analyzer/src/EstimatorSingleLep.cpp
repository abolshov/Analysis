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

ArrF_t<ESTIM_OUT_SZ> EstimatorSingleLep::EstimateCombination(VecLVF_t const& particles, std::vector<Float_t> const& jet_res, ULong64_t evt_id, JetComb const& comb)
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
            m_res_mass->Fill(mass, weight);
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
        return res;
    }
    
    return res;
}

OptArrF_t<ESTIM_OUT_SZ> EstimatorSingleLep::EstimateMass(VecLVF_t const& jets, std::vector<Float_t> const& resolutions, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id)
{
    VecLVF_t particles(static_cast<size_t>(ObjSL::count));
    particles[static_cast<size_t>(ObjSL::lep)] = leptons[static_cast<size_t>(Lep::lep1)];
    particles[static_cast<size_t>(ObjSL::met)] = met;

    std::vector<ArrF_t<ESTIM_OUT_SZ>> results;
    std::vector<Float_t> masses, inv_widths, integrals, peaks;
    size_t num_bjets = jets.size() < NUM_BEST_BTAG ? jets.size() : NUM_BEST_BTAG;

    // in case when aggregating estimations for all combinations into one distribution
    // name contains only event ID
    if (m_aggr_mode == AggregationMode::Event)
    {
        m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu", evt_id));
    }

    std::unordered_set<size_t> used;
    for (size_t bj1_idx = 0; bj1_idx < num_bjets; ++bj1_idx)
    {
        used.insert(bj1_idx);
        for (size_t bj2_idx = bj1_idx + 1; bj2_idx < num_bjets; ++bj2_idx)
        {
            used.insert(bj2_idx);
            if (jets[bj1_idx].Pt() > jets[bj2_idx].Pt())
            {
                particles[static_cast<size_t>(ObjSL::bj1)] = jets[bj1_idx];
                particles[static_cast<size_t>(ObjSL::bj2)] = jets[bj2_idx];
            }
            else 
            {
                particles[static_cast<size_t>(ObjSL::bj1)] = jets[bj2_idx];
                particles[static_cast<size_t>(ObjSL::bj2)] = jets[bj1_idx];
            }

            for (size_t lj1_idx = 0; lj1_idx < jets.size(); ++lj1_idx)
            {
                if (used.count(lj1_idx))
                {
                    continue;
                }
                used.insert(lj1_idx);

                for (size_t lj2_idx = lj1_idx + 1; lj2_idx < jets.size(); ++lj2_idx)
                {
                    if (used.count(lj2_idx))
                    {
                        continue;
                    }
                    used.insert(lj2_idx);

                    if (jets[lj1_idx].Pt() > jets[lj2_idx].Pt())
                    {
                        particles[static_cast<size_t>(ObjSL::lj1)] = jets[lj1_idx];
                        particles[static_cast<size_t>(ObjSL::lj2)] = jets[lj2_idx];
                    }
                    else 
                    {
                        particles[static_cast<size_t>(ObjSL::lj1)] = jets[lj2_idx];
                        particles[static_cast<size_t>(ObjSL::lj2)] = jets[lj1_idx];
                    }

                    JetComb comb{};
                    comb.b1 = bj1_idx;
                    comb.b2 = bj2_idx;
                    comb.q1 = lj1_idx;
                    comb.q2 = lj2_idx;
                    #ifdef EXPERIMENTAL 
                        ArrF_t<ESTIM_OUT_SZ> comb_result = Experimental::SL::EstimateCombination(particles, 
                                                                                                 m_pdf_1d,
                                                                                                 m_pdf_2d,
                                                                                                 m_res_mass,
                                                                                                 m_prg,                                                                        
                                                                                                 evt_id, 
                                                                                                 comb.ToString(Channel::SL));
                    #else 
                        ArrF_t<ESTIM_OUT_SZ> comb_result = EstimateCombination(particles, resolutions, evt_id, comb);
                    #endif

                    if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
                    {
                        results.push_back(comb_result);
                        masses.push_back(comb_result[static_cast<size_t>(EstimOut::mass)]);
                        inv_widths.push_back(1.0/comb_result[static_cast<size_t>(EstimOut::width)]);
                        peaks.push_back(comb_result[static_cast<size_t>(EstimOut::peak_value)]);
                        integrals.push_back(comb_result[static_cast<size_t>(EstimOut::integral)]);
                    }

                    // reset m_res_mass and use "clean" hist to build distribution for each combination
                    // else keep filling histogram and only reset it when moving to another event
                    if (m_aggr_mode == AggregationMode::Combination)
                    {
                        ResetHist(m_res_mass);
                    }
                    used.erase(lj2_idx);
                }
                used.erase(lj1_idx);
            }
            used.erase(bj2_idx);
        }
        used.erase(bj1_idx);
    }

    if (!results.empty())
    {
        if (m_aggr_mode == AggregationMode::Combination)
        {
            // MinMaxTransform(integrals.begin(), integrals.end());
            // MinMaxTransform(inv_widths.begin(), inv_widths.end());
            // MinMaxTransform(peaks.begin(), peaks.end());
            // MinMaxTransform(masses.begin(), masses.end());

            // Float_t max_metric = 0.0f;
            // for (size_t c = 0; c < results.size(); ++c)
            // {
            //     // Float_t metric = (integrals[c] + inv_widths[c] + peaks[c] + masses[c])/4.0;
            //     Float_t metric = (integrals[c] + inv_widths[c] + peaks[c])/3.0;
            //     if (metric > max_metric)
            //     {
            //         max_metric = metric;
            //         choice = c;
            //     }
            // }

            auto it = std::max_element(masses.begin(), masses.end());
            size_t choice = it - masses.begin();

            // auto it = std::ranges::max_element(integrals);
            // auto it = std::max_element(integrals.begin(), integrals.end());
            // size_t choice = it - integrals.begin();

            // auto it = std::max_element(peaks.begin(), peaks.end());
            // size_t choice = it - peaks.begin();            

            ResetHist(m_res_mass);
            return std::make_optional<ArrF_t<ESTIM_OUT_SZ>>(results[choice]);
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