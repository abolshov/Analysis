#include "EstimatorTools.hpp"
#include "RealQuadEqn.hpp"

#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TVector2.h"

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
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);
    [[maybe_unused]] auto [res1, res2] = lj_pt_res;
    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double eta = rg.Uniform(-6, 6);

        [[maybe_unused]] double hadW_mass = (j1 + j2).M();
        [[maybe_unused]] double lepW_mass = rg.Uniform(0, mh - hadW_mass);

        double dpt_1 = rg.Gaus(0, res1);
        bool valid_corr_1 = j1.Pt() + dpt_1 > 0.0;
        double dpx_1 = valid_corr_1 ? j1.Px() : 0.0;
        double dpy_1 = valid_corr_1 ? j1.Py() : 0.0;
        if (valid_corr_1)
        {
            j1.SetPtEtaPhiM(j1.Pt() + dpt_1, j1.Eta(), j1.Phi(), j1.M());
            
            dpx_1 -= j1.Px();
            dpx_1 *= -1.0;
            
            dpy_1 -= j1.Py();
            dpy_1 *= -1.0;
        }

        double dpt_2 = rg.Gaus(0, res2);
        bool valid_corr_2 = j2.Pt() + dpt_2 > 0.0;
        double dpx_2 = valid_corr_2 ? j2.Px() : 0.0;
        double dpy_2 = valid_corr_2 ? j2.Py() : 0.0;
        if (valid_corr_2)
        {
            j2.SetPtEtaPhiM(j2.Pt() + dpt_2, j2.Eta(), j2.Phi(), j2.M());

            dpx_2 -= j2.Px();
            dpx_2 *= -1.0;
            
            dpy_2 -= j2.Py();
            dpy_2 *= -1.0;
        }

        OptionalPair b_jet_resc = JetRescFact(b1, b2, b_resc_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() >= b2.Pt());
            auto [c1, c2] = b_jet_resc.value();
            TLorentzVector bb1 = c1*b1;
            TLorentzVector bb2 = c2*b2;

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
