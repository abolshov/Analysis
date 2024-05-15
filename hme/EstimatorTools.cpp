#include "EstimatorTools.hpp"

#include "TString.h"
#include "TCanvas.h"

OptPair AnalyticalMass(std::vector<TLorentzVector> const& particles)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    TLorentzVector vis1(l);
    vis1 += b1;
    vis1 += b2;
    vis1 += j1;
    vis1 += j2;

    double a = HIGGS_MASS*HIGGS_MASS - vis1.M()*vis1.M() + 2.0*vis1.Px()*met.Px() + 2.0*vis1.Py()*met.Py();
    double A = 4.0*(vis1.E()*vis1.E() - vis1.Pz()*vis1.Pz());
    double B = -4.0*a*vis1.Pz();
    double C = -4.0*vis1.E()*vis1.E()*(met.Px()*met.Px() + met.Py()*met.Py()) - a*a;
    double delta = B*B - 4.0*A*C;

    if (delta < 0.0) 
    {
        return std::nullopt;
    }

    double pz_1 = (-B + std::sqrt(delta))/(2.0*A);
    double pz_2 = (-B - std::sqrt(delta))/(2.0*A);

    TLorentzVector nu1(met.Px(), met.Py(), pz_1, met.E());
    TLorentzVector nu2(met.Px(), met.Py(), pz_2, met.E());

    TLorentzVector vis2(vis1);
    vis1 += nu1;
    vis2 += nu2;

    return std::make_optional<std::pair<double, double>>(vis1.M(), vis2.M());
}

OptPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass)
{
    if (j1.Pt() < j2.Pt()) 
    {
        std::swap(j1, j2);
    }

    for (int i = 0; i < N_ATTEMPTS; ++i)
    {
        double c1 = pdf->GetRandom();
    
        double x1 = j2.M2();
        double x2 = 2*c1*(j1*j2);
        double x3 = c1*c1*j1.M2() - mass*mass;

        double discriminant = (x2*x2 - 4.0*x1*x3);

        if (x2 < 0.0)
        {
            continue;
        }
        if (discriminant < 0.0 || std::abs(x1) < TOL)
        {
            continue;
        }

        double c2 = (-1.0*x2 + sqrt(discriminant))/(2*x1);

        if (c2 < 0.0)
        {
            continue;
        }

        return std::make_optional<std::pair<double, double>>(c1, c2);
    }

    return std::nullopt;
}

OptTLorVec ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, double mh, double eta)
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

OptTLorVec ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, double mh, double eta, double phi)
{
    TLorentzVector vis(l);
    vis += j1;
    vis += j2;

    double pt = (mh*mh - vis*vis)/(2*(vis.E()*std::cosh(eta) - vis.Px()*std::cos(phi) - vis.Py()*std::sin(phi) - vis.Pz()*std::sinh(eta)));

    if (std::isinf(pt) || std::isnan(pt) || pt < 0)
    {
        return std::nullopt;
    }

    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, 0.0);
    return std::make_optional<TLorentzVector>(std::move(v));
}

OptPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    rg.SetSeed(0);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    double W_dijet_mass = (j1 + j2).M();
    // if (W_dijet_mass > mh)
    // {
    //     return std::nullopt;
    // }

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double eta = rg.Uniform(-6, 6);
        OptPair b_jet_resc = JetRescFact(b1, b2, b_resc_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() > b2.Pt());
            auto [c1, c2] = b_jet_resc.value();
            TLorentzVector bb1 = c1*b1;
            TLorentzVector bb2 = c2*b2;
            double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px();
            double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py();

            TLorentzVector met_corr(met_corr_px, met_corr_py, 0.0, 0.0);
            OptTLorVec nu = ComputeNu(l, j1, j2, met_corr, mh, eta);

            if (nu)
            {
                TLorentzVector tmp(l);
                tmp += nu.value();

                double W_lep_mass = tmp.M();
                if (W_lep_mass > W_dijet_mass)
                {
                    // here W_dijet_mass is off-shell
                    if (W_dijet_mass > mh/2)
                    {
                        ++failed_iter;
                        continue;
                    }
                }
                else 
                {
                    // here W_lep_mass is off-shell
                    if (W_lep_mass > mh/2)
                    {
                        ++failed_iter;
                        continue;
                    }
                }

                if (W_lep_mass + W_dijet_mass > mh)
                {
                    ++failed_iter;
                    continue;
                } 

                tmp += j1;
                tmp += j2;  
                
                double H_WW_mass = tmp.M();

                if (std::abs(H_WW_mass - mh) > 1.0)
                {
                    ++failed_iter;
                    continue;
                }

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
