#include "tools.hpp"
#include "Tracer.hpp"
#include "constants.hpp"
#include "Utils.hpp"

#include <iostream>
#include <memory>
#include <random>
#include <cassert>
#include <algorithm>
#include <cmath>

extern TLorentzVector const zero;

void save::save_1d_dist(TH1F* dist, 
                        std::string const& name,
                        std::string const& title)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

    dist->GetXaxis()->SetTitle(title.c_str());
    dist->SetLineWidth(3);
    dist->SetStats(1);
    dist->DrawNormalized();
    c1->SaveAs((name + ".png").c_str());

    delete c1;
}

void save::save_2d_dist(TH2F* dist, 
                        std::string const& name,
                        std::string const& title_x,
                        std::string const& title_y)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();
    TStyle* gStyle = new TStyle();
    gStyle->SetPalette(kRainBow);

    dist->GetXaxis()->SetTitle(title_x.c_str());
    dist->GetYaxis()->SetTitle(title_y.c_str());
    dist->SetLineWidth(3);
    dist->Draw("colz");
    c1->SaveAs((name + ".png").c_str());

    delete gStyle;
    delete c1;
}

void save::save_1d_stack(std::vector<TH1F*> const& distrs,
                         std::vector<std::string> const& legends,
                         std::string const& name,
                         std::string const& title,
                         std::string const& axis_label)
{
    if (distrs.size() != legends.size())
    {
        std::cout << "number of legends and histograms do not match!" << std::endl;
        return;
    }
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

    THStack* stack = new THStack("stack", title.c_str());
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (int i = 0; i < static_cast<int>(distrs.size()); ++i)
    {
        distrs[i]->SetLineWidth(3);
        int line_color = i + 1 < 5 ? i + 1 : i + 2;
        distrs[i]->SetLineColor(line_color);
        stack->Add(distrs[i]);
        legend->AddEntry(distrs[i], legends[i].c_str(), "l");
    }

    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle(axis_label.c_str());
    legend->Draw();
    c1->SaveAs((name + ".png").c_str());

    delete stack;
    delete legend;
    delete c1;
}

std::pair<double, double> save::save_fit(TH1F* dist, std::string const& name, std::string const& title)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();
    
    dist->SetStats(0);

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    dist->GetXaxis()->SetTitle(title.c_str());
    dist->SetLineWidth(3);
    dist->SetLineColor(1);
    dist->Draw();

    auto fitFunc = new TF1("fitFunc", "landau", 0.0, 2000.0);
    fitFunc->SetParameter(1, 400);
    fitFunc->SetParameter(2, 5);
    fitFunc->SetParameter(0, dist->GetEntries());
    dist->Fit(fitFunc, "R");

    int ndf = fitFunc->GetNDF();
    double chi2 = fitFunc->GetChisquare();
    std::cout << "chi2/ndf = " << chi2 << "/" << ndf << " = " << chi2/ndf << std::endl;

    legend->AddEntry(dist, "Histogram", "l");
    legend->AddEntry(fitFunc, "Fit", "l");
    legend->Draw();

    c1->SaveAs((name + ".png").c_str());

    auto result = std::make_pair(fitFunc->GetParameter(1), fitFunc->GetParameter(2));

    delete legend;
    delete fitFunc;
    delete c1;

    return result;
}

void jet::pt_order(TLorentzVector& p1, TLorentzVector& p2)
{
    if (p1.Pt() < p2.Pt())
    {
        TLorentzVector tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
}

bool jet::is_offshell(TLorentzVector const& p1,
                      TLorentzVector const& p2,
                      TLorentzVector const& l,
                      TLorentzVector const& nu)
{
    return ((p1 + p2).M() < (l + nu).M());
}

std::pair<float, float> jet::compute_resc_factors(TLorentzVector& p1, 
                                                  TLorentzVector& p2,
                                                  TH1F* lead_jet_pdf, 
                                                  TH1F* mass_pdf)
{
    jet::pt_order(p1, p2);
    
    float mass = mass_pdf->GetRandom();

    for (int i = 0; i < 1000; ++i)
    {
        float lead_rescale_factor = lead_jet_pdf->GetRandom();
        
        // mass = mass_pdf->GetRandom(); 

        float x1 = p2.M2();
        float x2 = 2*lead_rescale_factor*(p1*p2);
        float x3 = lead_rescale_factor*lead_rescale_factor*p1.M2() - mass*mass;

        float discriminant = (x2 * x2 - 4 * x1 * x3);

        if (x2 < 0.0)
        {
            continue;
        }
        if (discriminant < 0.0 or x1 == 0.0)
        {
            continue;
        }

        float trail_rescale_factor = (-x2 + sqrt(discriminant)) / (2 * x1);

        if (trail_rescale_factor < 0.0)
        {
            continue;
        }

        return std::make_pair(lead_rescale_factor, trail_rescale_factor);
    }

    return FAIL;
}

float hme::analytical(std::vector<TLorentzVector> const& particles)
{
    float mass = -1.0;

    float mh = 125.0;
    TLorentzVector b1 = particles[0];
    TLorentzVector b2 = particles[1];
    TLorentzVector j1 = particles[2];
    TLorentzVector j2 = particles[3];
    TLorentzVector l = particles[4];
    TLorentzVector met = particles[5];

    TLorentzVector vis = (l + j1 + j2 + b1 + b2);
    float a = mh*mh - vis.M()*vis.M() + 2.0*vis.Px()*met.Px() + 2.0*vis.Py()*met.Py();
    float A = 4.0*(vis.E()*vis.E() - vis.Pz()*vis.Pz());
    float B = -4.0*a*vis.Pz();
    float C = -4.0*vis.E()*vis.E()*(met.Px()*met.Px() + met.Py()*met.Py()) - a*a;
    float delta = B*B - 4.0*A*C;

    if (delta < 0.0f) return mass;

    float pz_1 = (-B + sqrt(delta))/(2.0*A);
    float pz_2 = (-B - sqrt(delta))/(2.0*A);

    TLorentzVector nu1, nu2;
    nu1.SetPxPyPzE(met.Px(), met.Py(), pz_1, met.E());
    nu2.SetPxPyPzE(met.Px(), met.Py(), pz_2, met.E());

    TLorentzVector h_nu1 = (l + nu1 + j1 + j2);
    TLorentzVector h_nu2 = (l + nu2 + j1 + j2);

    TLorentzVector nu;
    if (abs(h_nu1.M() - 125.0) < abs(h_nu2.M() - 125.0))
    {
        nu = nu1;
    }
    else
    {
        nu = nu2;
    }

    // std::cout << "From Simple:" << std::endl;
    // std::cout << "----corrected dijet mass = " << (j1 + j2).M() << std::endl;
    // std::cout << "----leptonic W mass = " << (l + nu).M() << std::endl;

    float tmp_hh_mass = (l + nu + j1 + j2 + b1 + b2).M();
    mass = tmp_hh_mass;
    return mass;
}

TLorentzVector hme::NuFromLeptonicW_v1(float nu_eta, float nu_phi, TLorentzVector const& l, float mw)
{
    float deta = nu_eta - l.Eta();
    float dphi = nu_phi - l.Phi();
    float pt = mw*mw/(2*l.Pt()*(cosh(deta) - cos(dphi)));

    bool invalid = std::isinf(pt) || std::isnan(pt);
    float tmp_pt = invalid ? 0.0f : pt;
    float tmp_eta = invalid ? 0.0f : nu_eta;
    float tmp_phi = invalid ? 0.0f : nu_phi;

    TLorentzVector res;
    res.SetPtEtaPhiM(tmp_pt, tmp_eta, tmp_phi, 0.0f);
    return res;
}

TLorentzVector hme::NuFromLeptonicW_v2(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, float mh, int control)
{
    TLorentzVector vis = l + j1+ j2;
    TLorentzVector res;

    float fin_pt, fin_eta, fin_phi;

    float vis_pt = vis.Pt();
    float vis_m = vis.M();

    float tmp_px = met.Px();
    float tmp_py = met.Py();

    TLorentzVector tmp_vec;
    tmp_vec.SetPxPyPzE(sqrt(vis_pt*vis_pt + vis_m*vis_m), 0.0f, vis.Pz(), vis.E());

    fin_pt = met.Pt();

    float chdeta = (mh*mh + 2*tmp_px*vis.Px() + tmp_py*vis.Py() - vis_m*vis_m)/(2*tmp_vec.Pt()*fin_pt);
    if (chdeta < 1.0f)
    {
        res.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    }
    else
    {
        fin_phi = met.Phi();
        float deta = acosh(chdeta);
        fin_eta = (control == 1) ? (tmp_vec.Eta() - deta) : (tmp_vec.Eta() + deta);
        if (abs(fin_eta) > 7.0) 
        {
            res.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        }
        res.SetPtEtaPhiM(fin_pt, fin_eta, fin_phi, 0.0);
    }
    return res;
}

float hme::rand_sampl(std::vector<TLorentzVector> const& particles, std::vector<TH1F*> const& pdfs, int nIter, TRandom3& rg, int nbins, hme::MODE mode, bool correct_light_jets)
{
    TLorentzVector b1(particles[hme::GEN_PART::b1]);
    TLorentzVector b2(particles[hme::GEN_PART::b2]);
    TLorentzVector j1(particles[hme::GEN_PART::j1]);
    TLorentzVector j2(particles[hme::GEN_PART::j2]);
    TLorentzVector l(particles[hme::GEN_PART::l]);
    TLorentzVector met(particles[hme::GEN_PART::met]);
    TLorentzVector nu(particles[hme::GEN_PART::nu]);

    TH1F* lead_bjet_pdf = pdfs[0];
    TH1F* lead_on = pdfs[1];
    TH1F* lead_off = pdfs[2];
    TH1F* onshell_w_from_qq = pdfs[3];
    TH1F* offshell_w_from_qq = pdfs[4];
    TH1F* h_mass = pdfs[5];
    TH1F* nu_eta_pdf = pdfs[6];

    // std::unique_ptr<TH1F> hh_mass = std::make_unique<TH1F>("evt_hh_mass", "evt_hh_mass", 3000, 0.0, 3000.0);
    std::unique_ptr<TH1F> hh_mass = std::make_unique<TH1F>("evt_hh_mass", "evt_hh_mass", nbins, 0.0, 3000.0);
    float mh = h_mass->GetRandom();

    float hadronic_w = (j1 + j2).M();
    if (hadronic_w > mh) return -1.0f;

    bool is_offshell = hadronic_w < 60.0f;  
    // if hadronic W is ofshell, leptonic is onshell and vice versa
    // float leptonic_w = is_offshell ? onshell_w_from_qq->GetRandom() : offshell_w_from_qq->GetRandom();  
    float offshell_mass = offshell_w_from_qq->GetRandom(); 
    float onshell_mass = onshell_w_from_qq->GetRandom();

    rg.SetSeed(0);
    for (int i = 0; i < nIter; ++i)
    {
        float nu_eta = nu_eta_pdf->GetRandom();

        std::pair<float, float> light_jet_resc{1.0f, 1.0f};
        if (correct_light_jets)
        {
            if (is_offshell)
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_off, offshell_w_from_qq);
            }
            else
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_on, onshell_w_from_qq);
            }
        }


        std::pair<float, float> b_jet_resc = jet::compute_resc_factors(b1, b2, lead_bjet_pdf, h_mass);

        // std::pair<float, float> fail(-1.0, -1.0);
        // if (light_jet_resc == fail || b_jet_resc == fail)
        // if (b_jet_resc == fail)
        // {
        //     ++fails_counter;
        //     continue;
        // }

        float dpx = rg.Gaus(0, MET_SIGMA);
        float dpy = rg.Gaus(0, MET_SIGMA);

        float c1 = b_jet_resc.first;
        float c2 = b_jet_resc.second;

        TLorentzVector met_corr;

        float c3, c4;
        float met_px_corr = 0.0;
        float met_py_corr = 0.0;
        
        c3 = light_jet_resc.first;
        c4 = light_jet_resc.second;

        met_px_corr -= (c1 - 1)*b1.Px();
        met_px_corr -= (c2 - 1)*b2.Px();
        met_px_corr -= (c3 - 1)*j1.Px();
        met_px_corr -= (c4 - 1)*j2.Px();

        met_py_corr -= (c1 - 1)*b1.Py();
        met_py_corr -= (c2 - 1)*b2.Py();
        met_py_corr -= (c3 - 1)*j1.Py();
        met_py_corr -= (c4 - 1)*j2.Py();

        float met_px = met.Px();
        float met_py = met.Py();
        met_px += dpx;
        met_px += met_px_corr;
        met_py += dpy;
        met_py += met_py_corr;
        met_corr.SetPxPyPzE(met_px, met_py, met.Pz(), met.E());

        if (mode == hme::MODE::WeightedV1)
        {
            float tmp_hh_mass;
            TLorentzVector tmp_hh_momentum(b1);
            tmp_hh_momentum *= c1;
            tmp_hh_momentum += c2*b2;
            tmp_hh_momentum += c3*j1;
            tmp_hh_momentum += c4*j2;
            tmp_hh_momentum += l;

            // TLorentzVector tmp_nu = NuFromLeptonicW_v1(nu_eta, met_corr.Phi(), l, leptonic_w);
            // tmp_hh_momentum += tmp_nu;
            // tmp_hh_mass = tmp_hh_momentum.M();
            // if (tmp_hh_mass < 2.0f*mh) continue;
            // hh_mass->Fill(tmp_hh_mass);

            TLorentzVector nu_off = hme::NuFromLeptonicW_v1(nu_eta, met_corr.Phi(), l, offshell_mass);
            TLorentzVector nu_on = hme::NuFromLeptonicW_v1(nu_eta, met_corr.Phi(), l, onshell_mass);

            if (nu_off == zero && nu_on == zero) continue;
            if (nu_off != zero && nu_on != zero)
            {
                tmp_hh_mass = (tmp_hh_momentum + nu_off).M();
                // if (tmp_hh_mass < 2.0f*mh) continue;
                hh_mass->Fill(tmp_hh_mass, 0.5);

                tmp_hh_mass = (tmp_hh_momentum + nu_on).M();
                // if (tmp_hh_mass < 2.0f*mh) continue;
                hh_mass->Fill(tmp_hh_mass, 0.5);
            }
            else
            {
                if (nu_off != zero)
                {
                    tmp_hh_mass = (tmp_hh_momentum + nu_off).M();
                    // if (tmp_hh_mass < 2.0f*mh) continue;
                    hh_mass->Fill(tmp_hh_mass);
                }
                else
                {
                    tmp_hh_mass = (tmp_hh_momentum + nu_on).M();
                    // if (tmp_hh_mass < 2.0f*mh) continue;
                    hh_mass->Fill(tmp_hh_mass);
                }
            }
        }
        
        if (mode == hme::MODE::Original)
        {
            TLorentzVector tmp_nu;
            tmp_nu.SetPtEtaPhiM(met_corr.Pt(), nu_eta, met_corr.Phi(), 0.0);
            float tmp_hh_mass;
            TLorentzVector tmp_hh_momentum(b1);
            tmp_hh_momentum *= c1;
            tmp_hh_momentum += c2*b2;
            tmp_hh_momentum += c3*j1;
            tmp_hh_momentum += c4*j2;
            tmp_hh_momentum += l;
            tmp_hh_momentum += tmp_nu;
            tmp_hh_mass = tmp_hh_momentum.M();
            // if (tmp_hh_mass < 2.0f*mh) continue;
            hh_mass->Fill(tmp_hh_mass);
        }

        if (mode == hme::MODE::WeightedV2)
        {
            float tmp_hh_mass;
            TLorentzVector tmp_hh_momentum(b1);
            tmp_hh_momentum *= c1;
            tmp_hh_momentum += c2*b2;
            tmp_hh_momentum += c3*j1;
            tmp_hh_momentum += c4*j2;
            tmp_hh_momentum += l;

            TLorentzVector nu0 = hme::NuFromLeptonicW_v2(l, j1, j2, met_corr, mh, 0);
            // std::cout << "nu0 = ";
            // Print(nu0);
            TLorentzVector nu1 = hme::NuFromLeptonicW_v2(l, j1, j2, met_corr, mh, 1);
            // std::cout << "nu1 = ";
            // Print(nu1);

            if (nu0 == zero && nu1 == zero) continue;
            if (nu0 != zero && nu1 != zero)
            {
                tmp_hh_mass = (tmp_hh_momentum + nu0).M();
                // if (tmp_hh_mass < 2.0f*mh) continue;
                hh_mass->Fill(tmp_hh_mass, 0.5);
                // std::cout << "\tFilling tmp_hh_mass = " << tmp_hh_mass << "\n";

                tmp_hh_mass = (tmp_hh_momentum + nu1).M();
                // if (tmp_hh_mass < 2.0f*mh) continue;
                hh_mass->Fill(tmp_hh_mass, 0.5);
                // std::cout << "\tFilling tmp_hh_mass = " << tmp_hh_mass << "\n";
            }
            else
            {
                if (nu0 != zero)
                {
                    tmp_hh_mass = (tmp_hh_momentum + nu0).M();
                    // if (tmp_hh_mass < 2.0f*mh) continue;
                    hh_mass->Fill(tmp_hh_mass);
                    // std::cout << "\tFilling tmp_hh_mass = " << tmp_hh_mass << "\n";
                }
                else
                {
                    tmp_hh_mass = (tmp_hh_momentum + nu1).M();
                    // if (tmp_hh_mass < 2.0f*mh) continue;
                    hh_mass->Fill(tmp_hh_mass);
                    // std::cout << "\tFilling tmp_hh_mass = " << tmp_hh_mass << "\n";
                }
            }
        }

        // std::cout << "True W mass = " << (l + nu).M() << "\n";
        // std::cout << "nu = ";
        // Print(nu);
        // std::cout << "l = ";
        // Print(l);
        // std::cout << "met = ";
        // Print(met);
        // std::cout << "lmet mass = " << (l + met).M() << "\n";
        // std::cout << "tmp_nu = ";
        // Print(tmp_nu);
        // std::cout << "W mass = " << (l + tmp_nu).M() << "\n";

        // tmp_nu = NuFromLeptonicW_v1(nu_eta, met_corr.Phi(), l, offshell_mass);
        // std::cout << "offshell_mass = " << offshell_mass << "\n";
        // std::cout << "NuFromLeptonicW_v1 = ";
        // Print(tmp_nu);
        // std::cout << "W mass v1 = " << (l + tmp_nu).M() << "\n";

        // tmp_nu = NuFromLeptonicW_v1(nu_eta, met_corr.Phi(), l, onshell_mass);
        // std::cout << "onshell_mass = " << onshell_mass << "\n";
        // std::cout << "NuFromLeptonicW_v1 = ";
        // Print(tmp_nu);
        // std::cout << "W mass v1 = " << (l + tmp_nu).M() << "\n";

        // tmp_nu = NuFromLeptonicW_v2(l, j1, j2, met_corr, mh, 0);
        // std::cout << "NuFromLeptonicW_v2_1 = ";
        // Print(tmp_nu);
        // std::cout << "W mass v2_1 = " << (l + tmp_nu).M() << "\n";

        // tmp_nu = NuFromLeptonicW_v2(l, j1, j2, met_corr, mh, 1);
        // std::cout << "NuFromLeptonicW_v2_2 = ";
        // Print(tmp_nu);
        // std::cout << "W mass v2_1 = " << (l + tmp_nu).M() << "\n";

        // std::cout << "-------------------------------------------------\n";
    }

    if (hh_mass->GetEntries() == 0) return -1.0f;

    int binmax = hh_mass->GetMaximumBin(); 
    float evt_hh_mass = hh_mass->GetXaxis()->GetBinCenter(binmax);

    return evt_hh_mass;
}

float hme::anal_sampl(std::vector<TLorentzVector> const& particles, std::vector<TH1F*> const& pdfs, int nIter, TRandom3& rg, int nbins, bool correct_light_jets)
{
    TLorentzVector b1(particles[hme::GEN_PART::b1]);
    TLorentzVector b2(particles[hme::GEN_PART::b2]);
    TLorentzVector j1(particles[hme::GEN_PART::j1]);
    TLorentzVector j2(particles[hme::GEN_PART::j2]);
    TLorentzVector l(particles[hme::GEN_PART::l]);
    TLorentzVector met(particles[hme::GEN_PART::met]);
    TLorentzVector nu(particles[hme::GEN_PART::nu]);

    TH1F* lead_bjet_pdf = pdfs[0];
    TH1F* lead_on = pdfs[1];
    TH1F* lead_off = pdfs[2];
    TH1F* onshell_w_from_qq = pdfs[3];
    TH1F* offshell_w_from_qq = pdfs[4];
    TH1F* h_mass = pdfs[5];

    std::unique_ptr<TH1F> hh_mass = std::make_unique<TH1F>("evt_hh_mass", "evt_hh_mass", nbins, 0.0, 3000.0);

    float hadronic_w = (j1 + j2).M();
    float mh = h_mass->GetRandom();
    if (hadronic_w > mh) return -1.0f;
    bool is_offshell = hadronic_w < 60.0f;

    // sampling loop
    // sample b jet rescaling factors, correct jet momenta and met and pass corrected values to analytical
    for (int i = 0; i < nIter; ++i)
    {
        std::pair<float, float> light_jet_resc{1.0f, 1.0f};
        if (correct_light_jets)
        {
            if (is_offshell)
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_off, offshell_w_from_qq);
            }
            else
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_on, onshell_w_from_qq);
            }
        }

        std::pair<float, float> b_jet_resc = jet::compute_resc_factors(b1, b2, lead_bjet_pdf, h_mass);

        float dpx = rg.Gaus(0, MET_SIGMA);
        float dpy = rg.Gaus(0, MET_SIGMA);

        auto [c1, c2] = b_jet_resc;
        auto [c3, c4] = light_jet_resc;

        TLorentzVector met_corr;

        float met_px_corr = 0.0;
        float met_py_corr = 0.0;

        met_px_corr -= (c1 - 1)*b1.Px();
        met_px_corr -= (c2 - 1)*b2.Px();
        met_px_corr -= (c3 - 1)*j1.Px();
        met_px_corr -= (c4 - 1)*j2.Px();

        met_py_corr -= (c1 - 1)*b1.Py();
        met_py_corr -= (c2 - 1)*b2.Py();
        met_py_corr -= (c3 - 1)*j1.Py();
        met_py_corr -= (c4 - 1)*j2.Py();

        float met_px = met.Px();
        float met_py = met.Py();
        met_px += dpx;
        met_px += met_px_corr;
        met_py += dpy;
        met_py += met_py_corr;
        met_corr.SetPxPyPzE(met_px, met_py, met.Pz(), met.E());

        TLorentzVector bb1(c1*b1);
        TLorentzVector bb2(c2*b2);
        TLorentzVector jj1(c3*j1);
        TLorentzVector jj2(c4*j2);

        std::vector<TLorentzVector> parts = {bb1, bb2, jj1, jj2, l, met_corr};
        float tmp_mass = hme::analytical(parts);
        if (tmp_mass > 0.0)
        {
            hh_mass->Fill(tmp_mass);
        }
        else
        {
            continue;
        }
    }

    if (hh_mass->GetEntries() != 0)
    {
        int binmax = hh_mass->GetMaximumBin(); 
        float evt_hh_mass = hh_mass->GetXaxis()->GetBinCenter(binmax);
        return evt_hh_mass;
    }
    return -1.0f;
}