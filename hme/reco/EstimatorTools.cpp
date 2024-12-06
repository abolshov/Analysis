#include "EstimatorTools.hpp"
#include "RealQuadEqn.hpp"

#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TVector2.h"

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
                          std::unique_ptr<TH2F>& pdf_b1b2, 
                          std::vector<std::unique_ptr<TH1F>>& pdfs_1d,
                          TRandom3& rg, int evt, 
                          std::pair<double, double> lj_pt_res)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector lj1 = particles[PhysObj::wj1];
    TLorentzVector lj2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    std::unique_ptr<TH1F>& pdf_numet_pt = pdfs_1d[PDF::numet_pt];
    std::unique_ptr<TH1F>& pdf_numet_dphi = pdfs_1d[PDF::numet_dphi];
    std::unique_ptr<TH1F>& pdf_nulep_deta = pdfs_1d[PDF::nulep_deta];
    std::unique_ptr<TH1F>& pdf_hh_dphi = pdfs_1d[PDF::hh_dphi];
    std::unique_ptr<TH1F>& pdf_mbb = pdfs_1d[PDF::mbb];

    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);
    [[maybe_unused]] auto [res1, res2] = lj_pt_res;

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);

    [[maybe_unused]] int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        TLorentzVector j1 = GenerateResCorrection(lj1, rg, res1);
        TLorentzVector j2 = GenerateResCorrection(lj2, rg, res2);

        double dpx_1 = j1.Px() - lj1.Px();
        double dpx_2 = j2.Px() - lj2.Px();
        double dpy_1 = j1.Py() - lj1.Py();
        double dpy_2 = j2.Py() - lj2.Py();

        double eta = l.Eta() + pdf_nulep_deta->GetRandom(&rg);
        double dphi = pdf_numet_dphi->GetRandom(&rg);
        double met_fraction = pdf_numet_pt->GetRandom(&rg);

        double c1 = 1.0;
        double c2 = 1.0;
        pdf_b1b2->GetRandom2(c1, c2, &rg);

        int bin = pdf_b1b2->FindBin(c1, c2);
        double w = pdf_b1b2->GetBinContent(bin);

        TLorentzVector bb1 = b1;
        TLorentzVector bb2 = b2;

        bb1 *= c1;
        bb2 *= c2;

        TLorentzVector Hbb = bb1 + bb2;
        bin = pdf_mbb->FindBin(Hbb.M());
        w *= pdf_mbb->GetBinContent(bin);

        double smear_dpx = rg.Gaus(0.0, MET_SIGMA);
        double smear_dpy = rg.Gaus(0.0, MET_SIGMA);

        double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - dpx_1 - dpx_2 + smear_dpx;
        double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - dpy_1 - dpy_2 + smear_dpy;

        TLorentzVector met_corr(met_corr_px, met_corr_py, 0.0, 0.0);
        double pt = met_fraction*met_corr.Pt();
        double phi = met_corr.Phi() + dphi;

        TLorentzVector nu;
        nu.SetPtEtaPhiM(pt, eta, phi, 0.0);

        TLorentzVector Hww(l);
        Hww += nu;
        Hww += j1;
        Hww += j2;

        bin = pdf_mbb->FindBin(Hww.M());
        w *= pdf_mbb->GetBinContent(bin);

        double hh_dphi = Hww.DeltaPhi(Hbb);
        bin = pdf_hh_dphi->FindBin(hh_dphi);
        w *= pdf_hh_dphi->GetBinContent(bin);

        double X_mass = (Hww + Hbb).M();

        if (X_mass < 2*mh)
        {
            ++failed_iter;
            continue;
        }

        res_mass->Fill(X_mass, w);
    }

    if (res_mass->GetEntries() && res_mass->Integral() > 0.0)
    {
        int binmax = res_mass->GetMaximumBin(); 
        double estimated_mass = res_mass->GetXaxis()->GetBinCenter(binmax);
        // double success_rate = 1.0 - static_cast<double>(failed_iter)/N_ITER;
        // return std::make_optional<std::pair<double, double>>(estimated_mass, success_rate);
        return std::make_optional<std::pair<double, double>>(estimated_mass, res_mass->Integral());
    }

    return std::nullopt;
}