#include "EstimatorTools.hpp"
#include "RealQuadEqn.hpp"

#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TVector2.h"

#ifdef DEBUG
#include <sstream>
#include <fstream>
#endif

#ifdef PLOT
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#endif

TLorentzVector GenerateResCorrection(TLorentzVector const& v, TRandom3& rg, double res)
{
    TLorentzVector result;
    double dpt = rg.Gaus(0, res);
    double pt = v.Pt();
    while (pt + dpt < 0.0)
    {
        dpt = rg.Gaus(0, res);
    }
    result.SetPtEtaPhiM(pt + dpt, v.Eta(), v.Phi(), v.M());
    return result;
}

OptionalPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass)
{
    if (j1.Pt() <= j2.Pt()) 
    {
        std::swap(j1, j2);
    }

    for (int i = 0; i < N_ATTEMPTS; ++i)
    {
        double c1 = pdf->GetRandom();
    
        double x1 = j2.M2();
        double x2 = 2*c1*(j1*j2);
        double x3 = c1*c1*j1.M2() - mass*mass;

        std::vector<double> solutions = QuadEqn<double>(x1, x2, x3).Solutions();

        double c2 = 0;
        if (solutions.empty())
        {
            continue;
        }
        else if (solutions.size() == 1)
        {
            c2 = solutions[0];
        }
        else if (solutions.size() == 2)
        {
            c2 = solutions[0] > 0.0 ? solutions[0] : solutions[1];
        }

        if (c2 <= 0.0)
        {
            continue;
        }
        
        return std::make_optional<std::pair<double, double>>(c1, c2);
    }

    return std::nullopt;
}

std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, double mh, double eta)
{
    TLorentzVector vis(l);
    vis += j1;
    vis += j2;

    double phi = met.Phi();
    //m_h^2 = (j + j + l + nu)^2
    double pt = (mh*mh - vis*vis)/(2*(vis.E()*std::cosh(eta) - vis.Px()*std::cos(phi) - vis.Py()*std::sin(phi) - vis.Pz()*std::sinh(eta)));

    if (std::isinf(pt) || std::isnan(pt) || pt < 0)
    {
        return std::nullopt;
    }

    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, 0.0);
    return std::make_optional<TLorentzVector>(std::move(v));
}

std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& met, double mw, double numet_dphi)
{
    // assume that numet_dphi is in [-max_numet_dphi, max_numet_dphi]
    // numet_dphi is generated from a PDF
    // compute possible angle of nu knowing phi of met and possible correction (generated randomly)
    double nu_phi = TVector2::Phi_mpi_pi(met.Phi() + numet_dphi);
    double lep_phi = l.Phi();

    // angle between lepton and neutrino
    double dphi = TVector2::Phi_mpi_pi(lep_phi - nu_phi);

    double ch_deta = std::cos(dphi) + (mw*mw)/(2.0*met.Pt()*l.Pt());
    if (std::abs(ch_deta) < 1.0)
    {
        return std::nullopt;
    }

    double deta = std::acosh(ch_deta);
    double nu_eta = l.Eta() + deta;

    TLorentzVector res;
    res.SetPtEtaPhiM(met.Pt(), nu_eta, nu_phi, 0.0);
    return std::make_optional<TLorentzVector>(res);
}

OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt, std::pair<double, double> lj_pt_res)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector lj1 = particles[PhysObj::wj1];
    TLorentzVector lj2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);
    [[maybe_unused]] auto [res1, res2] = lj_pt_res;
    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    [[maybe_unused]] double hadW_mass = (lj1 + lj2).M();
    // [[maybe_unused]] double lepW_mass = rg.Uniform(0, mh - hadW_mass);

    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        TLorentzVector j1 = GenerateResCorrection(lj1, rg, res1);
        TLorentzVector j2 = GenerateResCorrection(lj2, rg, res2);

        double dpx_1 = j1.Px() - lj1.Px();
        double dpx_2 = j2.Px() - lj2.Px();
        double dpy_1 = j1.Py() - lj1.Py();
        double dpy_2 = j2.Py() - lj2.Py();

        double eta = rg.Uniform(-6, 6);

        // double dpt_1 = rg.Gaus(0, res1);
        // bool valid_corr_1 = j1.Pt() + dpt_1 > 0.0;
        // double dpx_1 = valid_corr_1 ? j1.Px() : 0.0;
        // double dpy_1 = valid_corr_1 ? j1.Py() : 0.0;
        // if (valid_corr_1)
        // {
        //     j1.SetPtEtaPhiM(j1.Pt() + dpt_1, j1.Eta(), j1.Phi(), j1.M());
            
        //     dpx_1 -= j1.Px();
        //     dpx_1 *= -1.0;
            
        //     dpy_1 -= j1.Py();
        //     dpy_1 *= -1.0;
        // }

        // double dpt_2 = rg.Gaus(0, res2);
        // bool valid_corr_2 = j2.Pt() + dpt_2 > 0.0;
        // double dpx_2 = valid_corr_2 ? j2.Px() : 0.0;
        // double dpy_2 = valid_corr_2 ? j2.Py() : 0.0;
        // if (valid_corr_2)
        // {
        //     j2.SetPtEtaPhiM(j2.Pt() + dpt_2, j2.Eta(), j2.Phi(), j2.M());

        //     dpx_2 -= j2.Px();
        //     dpx_2 *= -1.0;
            
        //     dpy_2 -= j2.Py();
        //     dpy_2 *= -1.0;
        // }

        // double mjj_corr = (j1 + j2).M();
        // if (mjj_corr > mh)
        // {
        //     ++failed_iter;
        //     continue;
        // }

        OptionalPair b_jet_resc = JetRescFact(b1, b2, b_resc_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() >= b2.Pt());
            auto [c1, c2] = b_jet_resc.value();
            TLorentzVector bb1 = b1;
            TLorentzVector bb2 = b2;

            bb1 *= c1;
            bb2 *= c2;

            double smear_dpx = rg.Gaus(0.0, MET_SIGMA);
            double smear_dpy = rg.Gaus(0.0, MET_SIGMA);

            double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - dpx_1 - dpx_2 + smear_dpx;
            double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - dpy_1 - dpy_2 + smear_dpy;
            // double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() + smear_dpx;
            // double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() + smear_dpy;

            TLorentzVector met_corr(met_corr_px, met_corr_py, 0.0, 0.0);
            std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, met_corr, mh, eta);

            if (nu)
            {
                TLorentzVector tmp(l);
                tmp += nu.value();
                tmp += j1;
                tmp += j2;  
                tmp += bb1;
                tmp += bb2;

                double X_mass = tmp.M();

                res_mass->Fill(X_mass);
            }
            else
            {
                ++failed_iter;
                continue;
            }
        }
        else
        {
            ++failed_iter;
            continue;
        }
    }

    if (res_mass->GetEntries())
    {
        int binmax = res_mass->GetMaximumBin(); 
        double estimated_mass = res_mass->GetXaxis()->GetBinCenter(binmax);
        double success_rate = 1.0 - static_cast<double>(failed_iter)/N_ITER;
        return std::make_optional<std::pair<double, double>>(estimated_mass, success_rate);
    }

    return std::nullopt;
}

OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, 
                          std::vector<std::unique_ptr<TH2F>>& pdfs_2d,
                          std::vector<std::unique_ptr<TH1F>>& pdfs_1d,
                          TRandom3& rg, int evt, int comb_id,
                          std::pair<double, double> lj_pt_res)
{
    TLorentzVector const& b1 = particles[PhysObj::bj1];
    TLorentzVector const& b2 = particles[PhysObj::bj2];
    TLorentzVector const& lj1 = particles[PhysObj::wj1];
    TLorentzVector const& lj2 = particles[PhysObj::wj2];
    TLorentzVector const& l = particles[PhysObj::lep];
    TLorentzVector const& met = particles[PhysObj::met];

    [[maybe_unused]] std::unique_ptr<TH2F>& pdf_b1b2 = pdfs_2d[PDF2::b1b2];
    [[maybe_unused]] std::unique_ptr<TH2F>& pdf_hh_dEtadPhi = pdfs_2d[PDF2::hh_dEtadPhi];
    [[maybe_unused]] std::unique_ptr<TH2F>& pdf_hh_pt_e = pdfs_2d[PDF2::hh_pt_e];

    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_numet_pt = pdfs_1d[PDF1::numet_pt];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_numet_dphi = pdfs_1d[PDF1::numet_dphi];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_nulep_deta = pdfs_1d[PDF1::nulep_deta];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_hh_dphi = pdfs_1d[PDF1::hh_dphi];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_mbb = pdfs_1d[PDF1::mbb];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_hh_deta = pdfs_1d[PDF1::hh_deta];
    [[maybe_unused]] std::unique_ptr<TH1F>& pdf_mww = pdfs_1d[PDF1::mww];

    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);
    auto [res1, res2] = lj_pt_res;

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d, comb %d", evt, comb_id), N_BINS, 200.0, MAX_MASS);

    #ifdef DEBUG
    std::stringstream log;
    log << "Event " << evt << ", raw values:\n" 
        << "bj1=(" << b1.Pt() << ", " << b1.Eta() << ", " << b1.Phi() << ", " << b1.M() << ")\n"
        << "bj2=(" << b2.Pt() << ", " << b2.Eta() << ", " << b2.Phi() << ", " << b2.M() << ")\n"
        << "lj1=(" << lj1.Pt() << ", " << lj1.Eta() << ", " << lj1.Phi() << ", " << lj1.M() << "), res=" << res1 << "\n"
        << "lj2=(" << lj2.Pt() << ", " << lj2.Eta() << ", " << lj2.Phi() << ", " << lj2.M() << "), res=" << res2 << "\n"
        << "lep=(" << l.Pt() << ", " << l.Eta() << ", " << l.Phi() << ", " << l.M() << ")\n"
        << "met=(" << met.Pt() << ", " << met.Eta() << ", " << met.Phi() << ", " << met.M() << ")\n"
        << "mbb=" << (b1 + b2).M() << "\n"
        << "mWhad=" << (lj1 + lj2).M() << "\n\n";
    #endif

    [[maybe_unused]] int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        #ifdef DEBUG
        log << "Iter " << i + 1 << ":\n";
        #endif
        TLorentzVector j1 = GenerateResCorrection(lj1, rg, res1);
        TLorentzVector j2 = GenerateResCorrection(lj2, rg, res2);

        #ifdef DEBUG
        log << "\tj1=(" << j1.Pt() << ", " << j1.Eta() << ", " << j1.Phi() << ", " << j1.M() << ")\n"
            << "\tj2=(" << j2.Pt() << ", " << j2.Eta() << ", " << j2.Phi() << ", " << j2.M() << ")\n"
            << "\thadW_mass=" << (j1 + j2).M() << "\n";
        #endif

        double dpx_1 = j1.Px() - lj1.Px();
        double dpx_2 = j2.Px() - lj2.Px();
        double dpy_1 = j1.Py() - lj1.Py();
        double dpy_2 = j2.Py() - lj2.Py();

        #ifdef DEBUG
        log << "\tdpx_1=" << dpx_1 << "\n"
            << "\tdpx_2=" << dpx_2 << "\n"
            << "\tdpy_1=" << dpy_1 << "\n"
            << "\tdpy_2=" << dpy_2 << "\n";
        #endif

        double eta = l.Eta() + pdf_nulep_deta->GetRandom(&rg);
        double dphi = pdf_numet_dphi->GetRandom(&rg);
        double met_fraction = pdf_numet_pt->GetRandom(&rg);

        #ifdef DEBUG
        log << "\teta=" << eta << "\n"
            << "\tdphi=" << dphi << "\n"
            << "\tmet_fraction=" << met_fraction << "\n";
        #endif

        double c1 = 1.0;
        double c2 = 1.0;
        pdf_b1b2->GetRandom2(c1, c2, &rg);

        #ifdef DEBUG
        log << "\tc1=" << c1 << ", c2=" << c2 << "\n";
        #endif

        TLorentzVector bb1 = b1;
        TLorentzVector bb2 = b2;

        bb1 *= c1;
        bb2 *= c2;

        TLorentzVector Hbb = bb1 + bb2;
    
        int bin = pdf_mbb->FindBin(Hbb.M());
        double w = pdf_mbb->GetBinContent(bin);

        #ifdef DEBUG
        log << "\tbb1=(" << bb1.Pt() << ", " << bb1.Eta() << ", " << bb1.Phi() << ", " << bb1.M() << ")\n"
            << "\tbb2=(" << bb2.Pt() << ", " << bb2.Eta() << ", " << bb2.Phi() << ", " << bb2.M() << ")\n"
            << "\tHbb_mass=" << Hbb.M() << "\n"
            << "\tdw=" << pdf_mbb->GetBinContent(pdf_mbb->FindBin(Hbb.M())) << ", w=" << w << "\n";
        #endif

        double smear_dpx = rg.Gaus(0.0, MET_SIGMA);
        double smear_dpy = rg.Gaus(0.0, MET_SIGMA);

        #ifdef DEBUG
        log << "\tsmear_dpx=" << smear_dpx << ", smear_dpy=" << smear_dpy << "\n";
        #endif

        double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - dpx_1 - dpx_2 + smear_dpx;
        double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - dpy_1 - dpy_2 + smear_dpy;

        TLorentzVector met_corr;
        double met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        double met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        met_corr.SetPtEtaPhiM(met_corr_pt, 0.0, met_corr_phi, 0.0);

        double pt = met_fraction*met_corr.Pt();
        double phi = TVector2::Phi_mpi_pi(met_corr.Phi() + dphi);

        TLorentzVector nu;
        nu.SetPtEtaPhiM(pt, eta, phi, 0.0);

        #ifdef DEBUG
        log << "\tmet_corr=(" << met_corr.Pt() << ", " << met_corr.Eta() << ", " << met_corr.Phi() << ", " << met_corr.M() << ")\n" 
            << "\tnu=(" << nu.Pt() << ", " << nu.Eta() << ", " << nu.Phi() << ", " << nu.M() << ")\n"; 
        #endif

        TLorentzVector Hww(l);
        Hww += nu;
        Hww += j1;
        Hww += j2;

        bin = pdf_mbb->FindBin(Hww.M());
        w *= pdf_mbb->GetBinContent(bin);

        double hh_dphi = Hbb.DeltaPhi(Hww);
        // bin = pdf_hh_dphi->FindBin(hh_dphi);
        // w *= pdf_hh_dphi->GetBinContent(bin);

        double hh_deta = Hbb.Eta() - Hww.Eta();
        // bin = pdf_hh_deta->FindBin(hh_deta);
        // w *= pdf_hh_deta->GetBinContent(bin);

        bin = pdf_hh_dEtadPhi->FindBin(hh_deta, hh_dphi);
        w *= pdf_hh_dEtadPhi->GetBinContent(bin);

        double r_hbb = Hbb.Pt()/Hbb.E();
        double r_hww = Hww.Pt()/Hww.E();
        bin = pdf_hh_pt_e->FindBin(r_hbb, r_hww);
        w *= pdf_hh_pt_e->GetBinContent(bin);

        #ifdef DEBUG
        log << "\tHww_mass=" << Hww.M() << ", dw=" << pdf_mww->GetBinContent(pdf_mww->FindBin(Hww.M())) << ", w=" << w << "\n"
            << "\thh_dphi=" << hh_dphi << ", hh_deta=" << hh_deta << ", dw=" << pdf_hh_dEtadPhi->GetBinContent(pdf_hh_dEtadPhi->FindBin(hh_deta)) << "\n"
            << "\tr_hbb=" << r_hbb << ", r_hww=" << r_hww << ", dw=" << pdf_hh_pt_e->GetBinContent(pdf_hh_pt_e->FindBin(hh_deta)) << ", w=" << w << "\n";
        #endif

        double X_mass = (Hww + Hbb).M();

        #ifdef DEBUG
        log << "\tX_mass=" << X_mass << "\n";
        log << "===========================================================\n";
        #endif

        if (X_mass < 2*mh)
        {
            ++failed_iter;
            continue;
        }


        res_mass->Fill(X_mass, w);
    }

    double integral = res_mass->Integral();
    if (res_mass->GetEntries() && integral > 0.0)
    {
        int binmax = res_mass->GetMaximumBin(); 
        double estimated_mass = res_mass->GetXaxis()->GetBinCenter(binmax);
        [[maybe_unused]] double max_content = res_mass->GetBinContent(binmax);
        // [[maybe_unused]] double width = InterquantileRange(res_mass);

        #ifdef PLOT
        auto canv = std::make_unique<TCanvas>("canv", "canv");
        canv->SetGrid();
        canv->SetTickx();
        canv->SetTicky();

        gStyle->SetOptStat();
        gStyle->SetStatH(0.25);

        res_mass->SetLineWidth(2);
        res_mass->Draw("hist");

        gPad->Update();

        auto stat_box = static_cast<TPaveStats*>(res_mass->GetListOfFunctions()->FindObject("stats"));   
        stat_box->AddText(Form("Width = %.2f", width));
        stat_box->AddText(Form("PeakX = %.2f", estimated_mass));
        stat_box->AddText(Form("PeakY = %.2e", max_content));
        stat_box->AddText(Form("Integral = %.2e", integral));
        stat_box->SetTextSize(0.03);
        stat_box->DrawClone();

        canv->SaveAs(Form("event/mass/evt_%d_comb_%d.png", evt, comb_id));
        #endif

        return std::make_optional<std::pair<double, double>>(estimated_mass, integral);
    }

    #ifdef DEBUG
    std::ofstream file(Form("event/debug/evt_%d_comb_%d.txt", evt, comb_id));
    file << log.str();
    file.close();
    #endif

    return std::nullopt;
}