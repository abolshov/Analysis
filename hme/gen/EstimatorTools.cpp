#include "EstimatorTools.hpp"
#include "RealQuadEqn.hpp"

#include <stdexcept>

#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"

OptionalPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass)
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
        // else 
        // {
        //     throw std::logic_error("Invalid (> 2) number of solutions of Higgs mass constraint");
        // }

        if (c2 <= 0.0)
        {
            continue;
        }
        
        return std::make_optional<std::pair<double, double>>(c1, c2);

        // double discriminant = (x2*x2 - 4.0*x1*x3);

        // if (x2 < 0.0)
        // {
        //     continue;
        // }
        // if (discriminant < 0.0 || std::abs(x1) < TOL)
        // {
        //     continue;
        // }

        // double c2 = (-1.0*x2 + sqrt(discriminant))/(2*x1);

        // if (c2 < 0.0)
        // {
        //     continue;
        // }

        // return std::make_optional<std::pair<double, double>>(c1, c2);
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

std::optional<TLorentzVector> ComputeNu(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, double mh, double eta, double phi)
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

OptionalPair EstimateMass(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt)
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

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double eta = rg.Uniform(-6, 6);
        OptionalPair b_jet_resc = JetRescFact(b1, b2, b_resc_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() > b2.Pt());
            auto [c1, c2] = b_jet_resc.value();
            TLorentzVector bb1 = c1*b1;
            TLorentzVector bb2 = c2*b2;
            double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px();
            double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py();

            TLorentzVector met_corr(met_corr_px, met_corr_py, 0.0, 0.0);
            std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, met_corr, mh, eta);

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

OptionalPair Experimental::EstimateMassIdealEtaPhi(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& b_resc_pdf, TRandom3& rg, int evt, std::pair<double, double> const& dir)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    auto const& [higgs_eta, higgs_phi] = dir;

    rg.SetSeed(42);
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
        // double phi = rg.Uniform(-3.1415, 3.1415);
        OptionalPair b_jet_resc = JetRescFact(b1, b2, b_resc_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() > b2.Pt());
            auto [c1, c2] = b_jet_resc.value();
            TLorentzVector bb1 = c1*b1;
            TLorentzVector bb2 = c2*b2;
            double met_corr_px = met.Px() - (c1 - 1)*b1.Px() - (c2 - 1)*b2.Px();
            double met_corr_py = met.Py() - (c1 - 1)*b1.Py() - (c2 - 1)*b2.Py();

            TLorentzVector met_corr(met_corr_px, met_corr_py, 0.0, 0.0);
            std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, met_corr, mh, eta);
            // std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, met, mh, eta);
            // std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, mh, eta, phi);

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

                TLorentzVector H_bb;
                TLorentzVector sum = bb1 + bb2;
                H_bb.SetPtEtaPhiM(sum.Pt(), higgs_eta, higgs_phi, sum.M());
                if (std::abs(H_bb.M() - mh) > 1.0)
                {
                    ++failed_iter;
                    continue;
                }
                tmp += H_bb;

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

OptionalPair Experimental::EstimateMassIdealNu(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& lead_b_pdf, TRandom3& rg, int evt, TLorentzVector const& nu)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    rg.SetSeed(42);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    double W_dijet_mass = (j1 + j2).M();

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        OptionalPair b_jet_resc = JetRescFact(b1, b2, lead_b_pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() > b2.Pt());
            auto [c1, c2] = b_jet_resc.value();

            TLorentzVector bb1 = c1*b1;
            TLorentzVector bb2 = c2*b2;

            TLorentzVector tmp(l);
            tmp += nu;

            // double W_lep_mass = tmp.M();
            // if (W_lep_mass > W_dijet_mass)
            // {
            //     // here W_dijet_mass is off-shell
            //     if (W_dijet_mass > mh/2)
            //     {
            //         ++failed_iter;
            //         continue;
            //     }
            // }
            // else 
            // {
            //     // here W_lep_mass is off-shell
            //     if (W_lep_mass > mh/2)
            //     {
            //         ++failed_iter;
            //         continue;
            //     }
            // }

            // if (W_lep_mass + W_dijet_mass > mh)
            // {
            //     ++failed_iter;
            //     continue;
            // } 

            tmp += j1;
            tmp += j2;  
            
            // double H_WW_mass = tmp.M();

            // if (std::abs(H_WW_mass - mh) > 1.0)
            // {
            //     ++failed_iter;
            //     continue;
            // }

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

    if (res_mass->GetEntries())
    {
        int binmax = res_mass->GetMaximumBin(); 
        double estimated_mass = res_mass->GetXaxis()->GetBinCenter(binmax);
        double success_rate = 1.0 - static_cast<double>(failed_iter)/N_ITER;

        // if (eff_lead_resc->GetEntries() && eff_subl_resc->GetEntries())
        // {
        //     auto c1 = std::make_unique<TCanvas>("c1", "c1");
        //     c1->SetGrid();
        //     c1->SetTickx();
        //     c1->SetTicky();

        //     binmax = eff_lead_resc->GetMaximumBin(); 
        //     double eff_c1 = eff_lead_resc->GetXaxis()->GetBinCenter(binmax);

        //     binmax = eff_subl_resc->GetMaximumBin(); 
        //     double eff_c2 = eff_subl_resc->GetXaxis()->GetBinCenter(binmax);

        //     std::cout << "\tEffective: c1=" << eff_c1 << ", c2=" << eff_c2 << "\n";

        //     eff_lead_resc->Draw();
        //     c1->SaveAs(Form("iter_plots/eff_lead_resc_%d.png", evt));
        //     eff_subl_resc->Draw();
        //     c1->SaveAs(Form("iter_plots/eff_subl_resc_%d.png", evt));
        // }

        return std::make_optional<std::pair<double, double>>(estimated_mass, success_rate);
    }

    return std::nullopt;
}

OptionalPair Experimental::EstimateMassIdealNu2dPDF(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH2F>& bjet_2d_pdf, TRandom3& rg, int evt, TLorentzVector const& nu)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    rg.SetSeed(42);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    if (b1.Pt() < b2.Pt())
    {
        std::swap(b1, b2);
    }

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double c1 = 1.0;
        double c2 = 1.0;

        bjet_2d_pdf->GetRandom2(c1, c2, &rg);
        TLorentzVector bb1 = c1*b1;
        TLorentzVector bb2 = c2*b2;

        TLorentzVector H_bb = bb1 + bb2;
        double H_bb_mass = H_bb.M();
        if (std::abs(H_bb_mass - mh) > 1.0)
        {
            ++failed_iter;
            continue;
        }

        TLorentzVector tmp(l);
        tmp += nu;

        // double W_lep_mass = tmp.M();
        // if (W_lep_mass > W_dijet_mass)
        // {
        //     // here W_dijet_mass is off-shell
        //     if (W_dijet_mass > mh/2)
        //     {
        //         ++failed_iter;
        //         continue;
        //     }
        // }
        // else 
        // {
        //     // here W_lep_mass is off-shell
        //     if (W_lep_mass > mh/2)
        //     {
        //         ++failed_iter;
        //         continue;
        //     }
        // }

        // if (W_lep_mass + W_dijet_mass > mh)
        // {
        //     ++failed_iter;
        //     continue;
        // } 

        tmp += j1;
        tmp += j2;  
        
        // double H_WW_mass = tmp.M();

        // if (std::abs(H_WW_mass - mh) > 1.0)
        // {
        //     ++failed_iter;
        //     continue;
        // }

        tmp += bb1;
        tmp += bb2;

        double X_mass = tmp.M();
        res_mass->Fill(X_mass);
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

OptionalPair Experimental::EstimateMassNew(std::vector<TLorentzVector> const& particles, std::vector<std::unique_ptr<TH1F>> const& pdfs, TRandom3& rg, int evt, TLorentzVector const& nu)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    std::unique_ptr<TH1F> const& lead_b_pt_pdf = pdfs[PDF::lead_bjet_pt_pdf];
    std::unique_ptr<TH1F> const& subl_b_pt_pdf = pdfs[PDF::subl_bjet_pt_pdf];

    std::unique_ptr<TH1F> const& lead_b_E_pdf = pdfs[PDF::lead_bjet_E_pdf];
    std::unique_ptr<TH1F> const& subl_b_E_pdf = pdfs[PDF::subl_bjet_E_pdf];

    rg.SetSeed(42);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    if (b1.Pt() < b2.Pt())
    {
        std::swap(b1, b2);
    }

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double c1_pt = lead_b_pt_pdf->GetRandom();
        double c2_pt = subl_b_pt_pdf->GetRandom();

        TLorentzVector bb1 = c1_pt*b1;
        TLorentzVector bb2 = c2_pt*b2;

        double c1_E = lead_b_E_pdf->GetRandom();
        double c2_E = subl_b_E_pdf->GetRandom();

        bb1.SetT(c1_E*bb1.E());
        bb2.SetT(c2_E*bb2.E());

        TLorentzVector H_bb = bb1 + bb2;
        double H_bb_mass = H_bb.M();
        if (std::abs(H_bb_mass - mh) > 1.0)
        {
            ++failed_iter;
            continue;
        }

        TLorentzVector tmp(l);
        tmp += nu;

        // double W_lep_mass = tmp.M();
        // if (W_lep_mass > W_dijet_mass)
        // {
        //     // here W_dijet_mass is off-shell
        //     if (W_dijet_mass > mh/2)
        //     {
        //         ++failed_iter;
        //         continue;
        //     }
        // }
        // else 
        // {
        //     // here W_lep_mass is off-shell
        //     if (W_lep_mass > mh/2)
        //     {
        //         ++failed_iter;
        //         continue;
        //     }
        // }

        // if (W_lep_mass + W_dijet_mass > mh)
        // {
        //     ++failed_iter;
        //     continue;
        // } 

        tmp += j1;
        tmp += j2;  
        
        // double H_WW_mass = tmp.M();

        // if (std::abs(H_WW_mass - mh) > 1.0)
        // {
        //     ++failed_iter;
        //     continue;
        // }

        tmp += bb1;
        tmp += bb2;

        double X_mass = tmp.M();
        res_mass->Fill(X_mass);
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

OptionalPair Experimental::EstimateMassIdealHbb(std::vector<TLorentzVector> const& particles, TRandom3& rg, int evt)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];
    TLorentzVector j1 = particles[PhysObj::wj1];
    TLorentzVector j2 = particles[PhysObj::wj2];
    TLorentzVector l = particles[PhysObj::lep];
    TLorentzVector met = particles[PhysObj::met];

    rg.SetSeed(42);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        double eta = rg.Uniform(-6, 6);
        std::optional<TLorentzVector> nu = ComputeNu(l, j1, j2, met, mh, eta);

        if (!nu)
        {
            ++failed_iter;
            continue;
        }

        double mass = (b1 + b2 + l + j1 + j2 + nu.value()).M();
        res_mass->Fill(mass);
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

OptionalPair Experimental::EstimateMassIdealHWW(std::vector<TLorentzVector> const& particles, std::unique_ptr<TH1F>& pdf, TRandom3& rg, int evt, TLorentzVector const& HtoWW)
{
    TLorentzVector b1 = particles[PhysObj::bj1];
    TLorentzVector b2 = particles[PhysObj::bj2];

    rg.SetSeed(42);
    double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);

    auto res_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", evt), N_BINS, 0.0, MAX_MASS);
    
    int failed_iter = 0;
    for (int i = 0; i < N_ITER; ++i)
    {
        OptionalPair b_jet_resc = JetRescFact(b1, b2, pdf, mh);
        if (b_jet_resc)
        {
            assert(b1.Pt() > b2.Pt());
            auto [c1, c2] = b_jet_resc.value();

            TLorentzVector bb1 = b1;
            bb1 *= c1;
            TLorentzVector bb2 = b2;
            bb2 *= c2;

            TLorentzVector tmp(HtoWW);
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

    if (res_mass->GetEntries())
    {
        int binmax = res_mass->GetMaximumBin(); 
        double estimated_mass = res_mass->GetXaxis()->GetBinCenter(binmax);
        double success_rate = 1.0 - static_cast<double>(failed_iter)/N_ITER;
        return std::make_optional<std::pair<double, double>>(estimated_mass, success_rate);
    }

    return std::nullopt;
}