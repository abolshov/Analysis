#include "Estimator.hpp"

#include <algorithm>
#include <unordered_set>

#include "TVector2.h"
#include "Math/GenVector/VectorUtil.h" // DeltaPhi
using ROOT::Math::VectorUtil::DeltaPhi;

#ifdef DEBUG 
#include <sstream>
#include <fstream>
#include "Analyzer.hpp"
#endif

#ifdef PLOT 
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#endif

Estimator::Estimator()
:   m_pdf_1d()
,   m_pdf_2d()
,   m_prg(std::make_unique<TRandom3>(SEED))
,   m_res_mass(std::make_unique<TH1F>("none", "none", N_BINS, MIN_MASS, MAX_MASS))
{}

EstimatorSingLep_Run3::EstimatorSingLep_Run3(TString const& file_name)
:   m_pdf_1d(pdf1d_sl_names.size())
,   m_pdf_2d(pdf2d_sl_names.size())
,   m_prg(std::make_unique<TRandom3>(SEED))
,   m_res_mass(std::make_unique<TH1F>("none", "none", N_BINS, MIN_MASS, MAX_MASS))
{
    TFile* pf = TFile::Open(file_name);
    Get1dPDFs(pf, m_pdf_1d, Channel::SL);
    Get2dPDFs(pf, m_pdf_2d, Channel::SL);
    pf->Close();
}

EstimatorDoubleLep_Run2::EstimatorDoubleLep_Run2(TString const& pdf_file_name)
:   Estimator()
{
    m_pdf_1d.reserve(pdf1d_dl_names.size());
    m_pdf_2d.reserve(pdf2d_dl_names.size());

    TFile* pf = TFile::Open(pdf_file_name);
    Get1dPDFs(pf, m_pdf_1d, Channel::DL);
    Get2dPDFs(pf, m_pdf_2d, Channel::DL);
    pf->Close();
}

std::array<Float_t, OUTPUT_SIZE> EstimatorSingLep_Run3::EstimateCombination(VecLVF_t const& particles, 
                                                                            std::pair<Float_t, Float_t> lj_pt_res,
                                                                            ULong64_t evt, TString const& comb_id)
{
    std::array<Float_t, OUTPUT_SIZE> res = {-1.0};

    LorentzVectorF_t const& bj1 = particles[static_cast<size_t>(ObjSL::bj1)];
    LorentzVectorF_t const& bj2 = particles[static_cast<size_t>(ObjSL::bj2)];
    LorentzVectorF_t const& lj1 = particles[static_cast<size_t>(ObjSL::lj1)];
    LorentzVectorF_t const& lj2 = particles[static_cast<size_t>(ObjSL::lj2)];
    LorentzVectorF_t const& lep = particles[static_cast<size_t>(ObjSL::lep)];
    LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjSL::met)];

    UHist_t<TH1F>& pdf_mbb = m_pdf_1d[static_cast<size_t>(PDF1_sl::mbb)];
    UHist_t<TH1F>& pdf_numet_pt = m_pdf_1d[static_cast<size_t>(PDF1_sl::numet_pt)];
    UHist_t<TH1F>& pdf_numet_dphi = m_pdf_1d[static_cast<size_t>(PDF1_sl::numet_dphi)];
    UHist_t<TH1F>& pdf_nulep_deta = m_pdf_1d[static_cast<size_t>(PDF1_sl::nulep_deta)];
    UHist_t<TH1F>& pdf_hh_dphi = m_pdf_1d[static_cast<size_t>(PDF1_sl::hh_dphi)];
    UHist_t<TH1F>& pdf_mww = m_pdf_1d[static_cast<size_t>(PDF1_sl::mww)];

    UHist_t<TH2F>& pdf_b1b2 = m_pdf_2d[static_cast<size_t>(PDF2_sl::b1b2)];
    UHist_t<TH2F>& pdf_hh_dEtadPhi = m_pdf_2d[static_cast<size_t>(PDF2_sl::hh_dEtadPhi)];
    UHist_t<TH2F>& pdf_hh_pt_e = m_pdf_2d[static_cast<size_t>(PDF2_sl::hh_pt_e)];

    Float_t mh = m_prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);
    auto [res1, res2] = lj_pt_res;

    m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt, comb_id.Data()));

    #ifdef DEBUG
        std::stringstream log;
        log << Analyzer::gen_truth_buf.str() << "\n";
        log << "reco values:\n" 
            << "\tbj1=(" << bj1.Pt() << ", " << bj1.Eta() << ", " << bj1.Phi() << ", " << bj1.M() << ")\n"
            << "\tbj2=(" << bj2.Pt() << ", " << bj2.Eta() << ", " << bj2.Phi() << ", " << bj2.M() << ")\n"
            << "\tlj1=(" << lj1.Pt() << ", " << lj1.Eta() << ", " << lj1.Phi() << ", " << lj1.M() << "), res=" << res1 << "\n"
            << "\tlj2=(" << lj2.Pt() << ", " << lj2.Eta() << ", " << lj2.Phi() << ", " << lj2.M() << "), res=" << res2 << "\n"
            << "\tlep=(" << lep.Pt() << ", " << lep.Eta() << ", " << lep.Phi() << ", " << lep.M() << ")\n"
            << "\tmet=(" << met.Pt() << ", " << met.Eta() << ", " << met.Phi() << ", " << met.M() << ")\n"
            << "\tmbb=" << (bj1 + bj2).M() << "\n"
            << "\tmWhad=" << (lj1 + lj2).M() << "\n\n";
    #endif

    [[maybe_unused]] int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        #ifdef DEBUG
            log << "Iter " << i + 1 << ":\n";
        #endif

        LorentzVectorF_t j1 = SamplePNetResCorr(lj1, m_prg, res1);
        LorentzVectorF_t j2 = SamplePNetResCorr(lj2, m_prg, res2);

        #ifdef DEBUG
            log << "\tj1=(" << j1.Pt() << ", " << j1.Eta() << ", " << j1.Phi() << ", " << j1.M() << ")\n"
                << "\tj2=(" << j2.Pt() << ", " << j2.Eta() << ", " << j2.Phi() << ", " << j2.M() << ")\n"
                << "\thadW_mass=" << (j1 + j2).M() << "\n";
        #endif

        Float_t dpx_1 = j1.Px() - lj1.Px();
        Float_t dpx_2 = j2.Px() - lj2.Px();
        Float_t dpy_1 = j1.Py() - lj1.Py();
        Float_t dpy_2 = j2.Py() - lj2.Py();

        #ifdef DEBUG
            log << "\tdpx_1=" << dpx_1 << "\n"
                << "\tdpx_2=" << dpx_2 << "\n"
                << "\tdpy_1=" << dpy_1 << "\n"
                << "\tdpy_2=" << dpy_2 << "\n";
        #endif

        Float_t eta = lep.Eta() + pdf_nulep_deta->GetRandom(m_prg.get());
        Float_t dphi = pdf_numet_dphi->GetRandom(m_prg.get());
        Float_t met_fraction = pdf_numet_pt->GetRandom(m_prg.get());

        #ifdef DEBUG
            log << "\teta=" << eta << "\n"
                << "\tdphi=" << dphi << "\n"
                << "\tmet_fraction=" << met_fraction << "\n";
        #endif

        LorentzVectorF_t b1 = bj1;
        LorentzVectorF_t b2 = bj2;
        Double_t c1 = 1.0;
        Double_t c2 = 1.0;
        pdf_b1b2->GetRandom2(c1, c2, m_prg.get());

        #ifdef DEBUG
            log << "\tc1=" << c1 << ", c2=" << c2 << "\n";
        #endif

        b1 *= c1;
        b2 *= c2;

        LorentzVectorF_t Hbb = b1 + b2;
    
        int bin = pdf_mbb->FindBin(Hbb.M());
        Float_t weight = pdf_mbb->GetBinContent(bin);

        #ifdef DEBUG
            log << "\tbb1=(" << b1.Pt() << ", " << b1.Eta() << ", " << b1.Phi() << ", " << b1.M() << ")\n"
                << "\tbb2=(" << b2.Pt() << ", " << b2.Eta() << ", " << b2.Phi() << ", " << b2.M() << ")\n"
                << "\tHbb_mass=" << Hbb.M() << "\n"
                << "\tdw=" << pdf_mbb->GetBinContent(pdf_mbb->FindBin(Hbb.M())) << ", weight=" << weight << "\n";
        #endif

        Float_t smear_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t smear_dpy = m_prg->Gaus(0.0, MET_SIGMA);

        #ifdef DEBUG
            log << "\tsmear_dpx=" << smear_dpx << ", smear_dpy=" << smear_dpy << "\n";
        #endif

        Float_t met_corr_px = met.Px() - (c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px() - dpx_1 - dpx_2 + smear_dpx;
        Float_t met_corr_py = met.Py() - (c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py() - dpy_1 - dpy_2 + smear_dpy;

        Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        LorentzVectorF_t met_corr(met_corr_pt, 0.0, met_corr_phi, 0.0);

        Float_t pt = met_fraction*met_corr.Pt();
        Float_t phi = TVector2::Phi_mpi_pi(met_corr.Phi() + dphi);

        LorentzVectorF_t nu(pt, eta, phi, 0.0);

        #ifdef DEBUG
            log << "\tmet_corr=(" << met_corr.Pt() << ", " << met_corr.Eta() << ", " << met_corr.Phi() << ", " << met_corr.M() << ")\n" 
                << "\tnu=(" << nu.Pt() << ", " << nu.Eta() << ", " << nu.Phi() << ", " << nu.M() << ")\n"; 
        #endif

        LorentzVectorF_t Hww = lep;
        Hww += nu;
        Hww += j1;
        Hww += j2;

        bin = pdf_mbb->FindBin(Hww.M());
        weight *= pdf_mbb->GetBinContent(bin);

        Float_t hh_dphi = DeltaPhi(Hbb, Hww);
        Float_t hh_deta = Hbb.Eta() - Hww.Eta();

        bin = pdf_hh_dEtadPhi->FindBin(hh_deta, hh_dphi);
        weight *= pdf_hh_dEtadPhi->GetBinContent(bin);

        Float_t r_hbb = Hbb.Pt()/Hbb.E();
        Float_t r_hww = Hww.Pt()/Hww.E();
        bin = pdf_hh_pt_e->FindBin(r_hbb, r_hww);
        weight *= pdf_hh_pt_e->GetBinContent(bin);

        #ifdef DEBUG
            log << "\tHww_mass=" << Hww.M() << ", dw=" << pdf_mbb->GetBinContent(pdf_mbb->FindBin(Hww.M())) << ", weight=" << weight << "\n"
                << "\thh_dphi=" << hh_dphi << ", hh_deta=" << hh_deta << ", dw=" << pdf_hh_dEtadPhi->GetBinContent(pdf_hh_dEtadPhi->FindBin(hh_deta)) << "\n"
                << "\tr_hbb=" << r_hbb << ", r_hww=" << r_hww << ", dw=" << pdf_hh_pt_e->GetBinContent(pdf_hh_pt_e->FindBin(hh_deta)) << ", weight=" << weight << "\n";
        #endif

        Float_t X_mass = (Hww + Hbb).M();

        #ifdef DEBUG
            log << "\tX_mass=" << X_mass << "\n";
            log << "===========================================================\n";
        #endif

        if (X_mass < 2.0*mh)
        {
            ++failed_iter;
            continue;
        }

        m_res_mass->Fill(X_mass, weight);
    }

    Float_t integral = m_res_mass->Integral();
    if (m_res_mass->GetEntries() && integral > 0.0)
    {
        int binmax = m_res_mass->GetMaximumBin(); 
        res[static_cast<size_t>(Output::mass)] = m_res_mass->GetXaxis()->GetBinCenter(binmax);
        res[static_cast<size_t>(Output::peak_val)] = m_res_mass->GetBinContent(binmax);
        res[static_cast<size_t>(Output::width)] = ComputeWidth(m_res_mass, Q16, Q84);
        res[static_cast<size_t>(Output::integral)] = integral;

        #ifdef DEBUG
            std::ofstream file(Form("event/debug/sl/evt_%llu_comb_%s.txt", evt, comb_id.Data()));
            file << log.str();
            file.close();
        #endif

        #ifdef PLOT
            auto canv = std::make_unique<TCanvas>("canv", "canv");
            canv->SetGrid();
            canv->SetTickx();
            canv->SetTicky();

            gStyle->SetOptStat();
            gStyle->SetStatH(0.25);

            m_res_mass->SetLineWidth(2);
            m_res_mass->Draw("hist");

            gPad->Update();

            auto stat_box = std::unique_ptr<TPaveStats>(static_cast<TPaveStats*>(m_res_mass->GetListOfFunctions()->FindObject("stats")));   
            stat_box->AddText(Form("Width = %.2f", res[static_cast<size_t>(Output::width)]));
            stat_box->AddText(Form("PeakX = %.2f", res[static_cast<size_t>(Output::mass)]));
            stat_box->AddText(Form("PeakY = %.2e", res[static_cast<size_t>(Output::peak_val)]));
            stat_box->AddText(Form("Integral = %.2e", res[static_cast<size_t>(Output::integral)]));
            stat_box->SetTextSize(0.03);
            stat_box->DrawClone();

            canv->SaveAs(Form("event/mass/sl/evt_%llu_comb_%s.png", evt, comb_id.Data()));
        #endif

        return res;
    }

    return res;
}

std::optional<Float_t> EstimatorSingLep_Run3::EstimateMass(VecLVF_t const& jets, 
                                                           VecLVF_t const& leptons, 
                                                           std::vector<Float_t> jet_resolutions, 
                                                           LorentzVectorF_t const& met, 
                                                           ULong64_t evt,
                                                           TString& chosen_comb)
{
    VecLVF_t particles(static_cast<size_t>(ObjSL::count));
    particles[static_cast<size_t>(ObjSL::lep)] = leptons[static_cast<size_t>(Lep::lep1)];
    particles[static_cast<size_t>(ObjSL::met)] = met;

    std::vector<Float_t> estimations;
    std::vector<Float_t> integrals;

    #ifdef DEBUG
        std::vector<TString> labels;
    #endif

    std::unordered_set<size_t> used;
    for (size_t bj1_idx = 0; bj1_idx < NUM_BEST_BTAG; ++bj1_idx)
    {
        used.insert(bj1_idx);
        for (size_t bj2_idx = bj1_idx + 1; bj2_idx < NUM_BEST_BTAG; ++bj2_idx)
        {
            used.insert(bj2_idx);

            // order jets such that first b jet has bigger pt and save their p4
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

                    particles[static_cast<size_t>(ObjSL::lj1)] = jets[lj1_idx];
                    particles[static_cast<size_t>(ObjSL::lj2)] = jets[lj2_idx];

                    std::pair<Float_t, Float_t> lj_res = {jet_resolutions[lj1_idx], jet_resolutions[lj2_idx]};

                    TString comb_label = Form("b%zub%zuq%zuq%zu", bj1_idx, bj2_idx, lj1_idx, lj2_idx);
                    auto comb_result = EstimateCombination(particles, lj_res, evt, comb_label);

                    // success: mass > 0
                    if (comb_result[static_cast<size_t>(Output::mass)] > 0.0)
                    {
                        estimations.push_back(comb_result[static_cast<size_t>(Output::mass)]);
                        integrals.push_back(comb_result[static_cast<size_t>(Output::integral)]);
                        
                        #ifdef DEBUG
                            labels.push_back(comb_label);
                        #endif
                    }

                    // clear the histogram to be reused 
                    Reset(m_res_mass);

                    // remove indices of light jets from used
                    used.erase(lj2_idx);
                }
                used.erase(lj1_idx);
            }
            used.erase(bj2_idx);
        }
        used.erase(bj1_idx);
    }

    // success: at least one combination produced an estimate of X->HH mass
    if (!estimations.empty())
    {
        auto it = std::max_element(integrals.begin(), integrals.end());
        int idx = it - integrals.begin();

        #ifdef DEBUG
            chosen_comb = labels[idx];
        #endif

        return std::make_optional<Float_t>(estimations[idx]);
    }

    return std::nullopt;
}

std::array<Float_t, OUTPUT_SIZE> EstimatorDoubleLep_Run2::EstimateCombination(VecLVF_t const& particles, 
                                                                              ULong64_t evt, 
                                                                              TString const& comb_id)
{
    std::array<Float_t, OUTPUT_SIZE> res = {-1.0};

    LorentzVectorF_t const& bj1 = particles[static_cast<size_t>(ObjDL::bj1)];
    LorentzVectorF_t const& bj2 = particles[static_cast<size_t>(ObjDL::bj2)];
    LorentzVectorF_t const& lep1 = particles[static_cast<size_t>(ObjDL::lep1)];
    LorentzVectorF_t const& lep2 = particles[static_cast<size_t>(ObjDL::lep2)];
    LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjDL::met)];

    UHist_t<TH1F>& pdf_b1 = m_pdf_1d[static_cast<size_t>(PDF1_dl::b1)];
    UHist_t<TH1F>& pdf_mw_onshell = m_pdf_1d[static_cast<size_t>(PDF1_dl::mw_onshell)];

    m_res_mass->SetNameTitle("X_mass", Form("X->HH mass: event %llu, comb %s", evt, comb_id.Data()));

    #ifdef DEBUG
        std::stringstream log;
        log << Analyzer::gen_truth_buf.str() << "\n";
        log << "reco values:\n" 
            << "\tbj1=(" << bj1.Pt() << ", " << bj1.Eta() << ", " << bj1.Phi() << ", " << bj1.M() << ")\n"
            << "\tbj2=(" << bj2.Pt() << ", " << bj2.Eta() << ", " << bj2.Phi() << ", " << bj2.M() << ")\n"
            << "\tlep1=(" << lep1.Pt() << ", " << lep1.Eta() << ", " << lep1.Phi() << ", " << lep1.M() << ")\n"
            << "\tlep2=(" << lep2.Pt() << ", " << lep2.Eta() << ", " << lep2.Phi() << ", " << lep2.M() << ")\n"
            << "\tmet=(" << met.Pt() << ", " << met.Eta() << ", " << met.Phi() << ", " << met.M() << ")\n"
            << "\tmbb=" << (bj1 + bj2).M() << "\n\n";
    #endif

    for (int i = 0; i < N_ITER; ++i)
    {
        #ifdef DEBUG
            log << "Iter " << i + 1 << ":\n";
        #endif

        Float_t eta_gen = m_prg->Uniform(-6, 6);
        Float_t phi_gen = m_prg->Uniform(-3.1415926, 3.1415926);
        Float_t mh = m_prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);
        Float_t mw = pdf_mw_onshell->GetRandom(m_prg.get());
        Float_t smear_dpx = m_prg->Gaus(0.0, MET_SIGMA);
        Float_t smear_dpy = m_prg->Gaus(0.0, MET_SIGMA);

        #ifdef DEBUG
            log << "\teta_gen=" << eta_gen << ", phi_gen=" << phi_gen << "\n"
                << "\tmh=" << mh << "mw=" << mw << "\n"
                << "\tsmear_dpx=" << smear_dpx << ", smear_dpy=" << smear_dpy << "\n";
        #endif

        auto [c1, c2] = ComputeJetResc(bj1, bj2, pdf_b1, mh);
        #ifdef DEBUG
            log << "\tc1=" << c1 << ", c2=" << c2 << "\n";
        #endif

        LorentzVectorF_t b1 = bj1;
        LorentzVectorF_t b2 = bj2;
        b1 *= c1;
        b2 *= c2;

        Float_t jet_resc_dpx = -(c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px();
        Float_t jet_resc_dpy = -(c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py();

        #ifdef DEBUG
            log << "\tjet_resc_dpx=" << jet_resc_dpx << ", jet_resc_dpy=" << jet_resc_dpy << "\n";
        #endif

        Float_t met_corr_px = met.Px() + jet_resc_dpx + smear_dpx;
        Float_t met_corr_py = met.Py() + jet_resc_dpy + smear_dpy;

        Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        LorentzVectorF_t met_corr(met_corr_pt, 0.0, met_corr_phi, 0.0);

        #ifdef DEBUG
            LogP4(log, met_corr, "met_corr");
            // log << "\tmet_corr=(" << met_corr.Pt() << ", " << met_corr.Eta() << ", " << met_corr.Phi() << ", " << met_corr.M() << ")\n";
        #endif

        LorentzVectorF_t Hbb = b1;
        Hbb += b2;

        #ifdef DEBUG
            LogP4(log, Hbb, "Hbb");
            log << "\tHbb_mass=" << Hbb.M() << "\n";
        #endif

        std::vector<Float_t> estimates;

        // two options: 
        // lep1 comes from onshell W and lep2 comes from offshell W
        // lep1 comes from offshell W and lep2 comes from onshell W
        // when neutrino is computed in each case again two options are possible: delta_eta is added or subtracted to nu
        // in total: 4 combinations; they are encoded in this for loop
        for (int j = 0; j < 4; ++j)
        {
            #ifdef DEBUG
                log << "\tsolution " << j + 1 << "/4\n";
            #endif

            LorentzVectorF_t l_offshell;
            LorentzVectorF_t l_onshell;
            
            int is_onshell = j/2;
            if (is_onshell == 0) 
            {
                l_onshell = lep1;
                l_offshell = lep2;
            } 
            else 
            {
                l_onshell = lep2;
                l_offshell = lep1;
            }

            auto nu_onshell = NuFromOnshellW(eta_gen, phi_gen, mw, l_onshell);  

            #ifdef DEBUG
                if (nu_onshell)
                {
                    LogP4(log, nu_onshell.value(), "nu_onshell");
                }
                else 
                {
                    log << "\tnu from onshell W not possible";
                }
            #endif

            if (!nu_onshell)
            {
                continue;
            }

            int is_offshell = j%2;
            auto nu_offshell = NuFromOffshellW(lep1, lep2, nu_onshell.value(), met_corr, is_offshell, mh);

            #ifdef DEBUG
                if (nu_offshell)
                {
                    LogP4(log, nu_offshell.value(), "nu_offshell");
                }
                else 
                {
                    log << "\tnu from offshell W not possible";
                }
            #endif

            if (!nu_offshell)
            {
                continue;
            }

            LorentzVectorF_t Hww = l_offshell;
            Hww += l_onshell;
            Hww += nu_onshell.value();
            Hww += nu_offshell.value(); 

            #ifdef DEBUG
                Float_t hh_dphi = DeltaPhi(Hww, Hbb);
                LogP4(log, Hww, "Hww");
                log << "\tHww_mass=" << Hww.M() << "\n"
                    << "\thh_dphi=" << hh_dphi << "\n";
            #endif

            if (std::abs(Hww.M() - mh) > 1.0)
            {
                continue;
            }

            Float_t mX = (Hbb + Hww).M();
            if (mX > 0.0)
            {
                estimates.push_back(mX);
            }
        }

        Float_t weight = estimates.empty() ? 0.0 : 1.0/estimates.size();
        for (auto est: estimates)
        {
            m_res_mass->Fill(est, weight);
        }

        #ifdef DEBUG
            log << "\n_solutions=" << estimates.size() << "\n"
                << "\tweight=" << weight << "\n";
        #endif
    }

    Float_t integral = m_res_mass->Integral();
    if (m_res_mass->GetEntries() && integral > 0.0)
    {
        int binmax = m_res_mass->GetMaximumBin(); 
        res[static_cast<size_t>(Output::mass)] = m_res_mass->GetXaxis()->GetBinCenter(binmax);
        res[static_cast<size_t>(Output::peak_val)] = m_res_mass->GetBinContent(binmax);
        res[static_cast<size_t>(Output::width)] = ComputeWidth(m_res_mass, Q16, Q84);
        res[static_cast<size_t>(Output::integral)] = integral;

        #ifdef DEBUG
            std::ofstream file(Form("event/debug/dl/evt_%llu_comb_%s.txt", evt, comb_id.Data()));
            file << log.str();
            file.close();
        #endif

        #ifdef PLOT
            auto canv = std::make_unique<TCanvas>("canv", "canv");
            canv->SetGrid();
            canv->SetTickx();
            canv->SetTicky();

            gStyle->SetOptStat();
            gStyle->SetStatH(0.25);

            m_res_mass->SetLineWidth(2);
            m_res_mass->Draw("hist");

            gPad->Update();

            auto stat_box = std::unique_ptr<TPaveStats>(static_cast<TPaveStats*>(m_res_mass->GetListOfFunctions()->FindObject("stats")));   
            stat_box->AddText(Form("Width = %.2f", res[static_cast<size_t>(Output::width)]));
            stat_box->AddText(Form("PeakX = %.2f", res[static_cast<size_t>(Output::mass)]));
            stat_box->AddText(Form("PeakY = %.2e", res[static_cast<size_t>(Output::peak_val)]));
            stat_box->AddText(Form("Integral = %.2e", res[static_cast<size_t>(Output::integral)]));
            stat_box->SetTextSize(0.03);
            stat_box->DrawClone();

            canv->SaveAs(Form("event/mass/dl/evt_%llu_comb_%s.png", evt, comb_id.Data()));
        #endif

        return res;
    }

    return res;
}

std::optional<Float_t> EstimatorDoubleLep_Run2::EstimateMass(VecLVF_t const& jets, 
                                                             VecLVF_t const& leptons, 
                                                             LorentzVectorF_t const& met, 
                                                             ULong64_t evt,
                                                             TString& chosen_comb)
{
    VecLVF_t particles(static_cast<size_t>(ObjDL::count));
    particles[static_cast<size_t>(ObjDL::lep1)] = leptons[static_cast<size_t>(Lep::lep1)];
    particles[static_cast<size_t>(ObjDL::lep2)] = leptons[static_cast<size_t>(Lep::lep2)];
    particles[static_cast<size_t>(ObjDL::met)] = met;

    std::vector<Float_t> estimations;
    std::vector<Float_t> integrals;

    #ifdef DEBUG
        std::vector<TString> labels;
    #endif

    for (size_t bj1_idx = 0; bj1_idx < NUM_BEST_BTAG; ++bj1_idx)
    {
        for (size_t bj2_idx = bj1_idx + 1; bj2_idx < NUM_BEST_BTAG; ++bj2_idx)
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
            auto comb_result = EstimateCombination(particles, evt, comb_label);

            // success: mass > 0
            if (comb_result[static_cast<size_t>(Output::mass)] > 0.0)
            {
                estimations.push_back(comb_result[static_cast<size_t>(Output::mass)]);
                integrals.push_back(comb_result[static_cast<size_t>(Output::integral)]);
                
                #ifdef DEBUG
                    labels.push_back(comb_label);
                #endif
            }

            // clear the histogram to be reused 
            Reset(m_res_mass);
        }
    }

    // success: at least one combination produced an estimate of X->HH mass
    if (!estimations.empty())
    {
        auto it = std::max_element(integrals.begin(), integrals.end());
        int idx = it - integrals.begin();

        #ifdef DEBUG
            chosen_comb = labels[idx];
        #endif

        return std::make_optional<Float_t>(estimations[idx]);
    }

    return std::nullopt;
}