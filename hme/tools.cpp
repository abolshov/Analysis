#include "tools.hpp"
#include "Tracer.hpp"

#include <iostream>
#include <memory>

bool random_hme_failed = false;
bool simpl_impr_hme_failed = false;
float random_hme_eff = 0.0;
float simpl_impr_hme_eff = 0.0;

void save::save_1d_dist(TH1F* dist, 
                        std::string const& name,
                        std::string const& title)
{
    TCanvas* c1 = new TCanvas("c1", "c1");

    dist->GetXaxis()->SetTitle(title.c_str());
    dist->SetLineWidth(3);
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

    THStack* stack = new THStack("stack", title.c_str());
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (size_t i = 0; i < distrs.size(); ++i)
    {
        distrs[i]->SetLineWidth(3);
        distrs[i]->SetLineColor(i + 1);
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

    for (size_t i = 0; i < 1000; ++i)
    {
        float lead_rescale_factor = lead_jet_pdf->GetRandom();
        
        // mass = mass_pdf->GetRandom(); 

        float x1 = p2.M2();
        float x2 = 2*lead_rescale_factor*(p1*p2);
        float x3 = lead_rescale_factor*lead_rescale_factor*p1.M2() - mass*mass;

        if (x2 < 0.0)
        {
            continue;
        }
        if ((x2 * x2 - 4 * x1 * x3) < 0.0 or x1 == 0.0)
        {
            continue;
        }

        float trail_rescale_factor = (-x2 + sqrt(x2 * x2 - 4 * x1 * x3)) / (2 * x1);

        if (trail_rescale_factor < 0.0)
        {
            continue;
        }

        return std::make_pair(lead_rescale_factor, trail_rescale_factor);
    }

    return std::make_pair(-1.0, -1.0);
}

float hme::hme_simplified(std::vector<TLorentzVector> const& particles)
{
    float mass = -1.0;

    static int fails_counter = 0;

    float mh = 125.0;
    TLorentzVector b1 = particles[0];
    TLorentzVector b2 = particles[1];
    TLorentzVector j1 = particles[2];
    TLorentzVector j2 = particles[3];
    TLorentzVector l = particles[4];
    TLorentzVector met = particles[5];

    // float universal_corr = mh/(b1+b2).M();

    // b1 *= universal_corr;
    // b2 *= universal_corr;

    TLorentzVector vis = (l + j1 + j2);
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

    float tmp_hh_mass = (l + nu + j1 + j2 + b1 + b2).M();
    mass = tmp_hh_mass;
    return mass;
}

float hme::hme_rand_sampl(std::vector<TLorentzVector> const& particles, std::vector<TH1F*> const& pdfs, int nIter, TRandom3& rg, bool light_on, bool uniform_eta)
{
    TLorentzVector b1 = particles[0];
    TLorentzVector b2 = particles[1];
    TLorentzVector j1 = particles[2];
    TLorentzVector j2 = particles[3];
    TLorentzVector l = particles[4];
    TLorentzVector nu = particles[5];
    TLorentzVector met = particles[6];

    TH1F* lead_bjet_pdf = pdfs[0];
    TH1F* lead_on = pdfs[1];
    TH1F* lead_off = pdfs[2];
    TH1F* onshell_w_from_qq = pdfs[3];
    TH1F* offshell_w_from_qq = pdfs[4];
    TH1F* h_mass = pdfs[5];
    TH1F* nu_eta_pdf = pdfs[6];

    float const met_sigma = 25.2;

    std::unique_ptr<TH1F> hh_mass = std::make_unique<TH1F>("evt_hh_mass", "evt_hh_mass", 120, 0.0, 3000.0);

    static int fails_counter = 0;

    rg.SetSeed(0);
    for (size_t i = 0; i < nIter; ++i)
    {
        float nu_eta;
        if (uniform_eta)
        {
            nu_eta = rg.Uniform(-6, 6);
        }
        else
        {
            nu_eta = nu_eta_pdf->GetRandom();
        }

        std::pair<float, float> light_jet_resc;
        if (light_on)
        {
            if (jet::is_offshell(j1, j2, l, nu))
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_off, offshell_w_from_qq);
            }
            else
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_on, onshell_w_from_qq);
            }
        }

        std::pair<float, float> b_jet_resc = jet::compute_resc_factors(b1, b2, lead_bjet_pdf, h_mass);

        std::pair<float, float> fail(-1.0, -1.0);

        // if (light_jet_resc == fail || b_jet_resc == fail)
        if (b_jet_resc == fail)
        {
            ++fails_counter;
            continue;
        }

        float dpx = rg.Gaus(0, met_sigma);
        float dpy = rg.Gaus(0, met_sigma);

        float c1 = b_jet_resc.first;
        float c2 = b_jet_resc.second;

        TLorentzVector met_corr;

        float c3, c4;
        float met_px_corr, met_py_corr;
        if (light_on)
        {
            c3 = light_jet_resc.first;
            c4 = light_jet_resc.second;
            met_px_corr = -(c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - (c3 - 1)*j1.Px() - (c4 - 1)*j2.Px();
            met_py_corr = -(c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - (c3 - 1)*j1.Py() - (c4 - 1)*j2.Py();
        }
        else
        {
            met_px_corr = -(c1 - 1)*b1.Px() - (c2 - 1)*b2.Px();
            met_py_corr = -(c1 - 1)*b1.Py() - (c2 - 1)*b2.Py();
        }

        met_corr.SetPxPyPzE(met.Px() + dpx + met_px_corr, met.Py() + dpy + met_py_corr, met.Pz(), met.E());

        TLorentzVector nu_corr;
        nu_corr.SetPtEtaPhiM(met_corr.Pt(), nu_eta, met_corr.Phi(), 0.0);

        float tmp_hh_mass;
        if (light_on)
        {
            tmp_hh_mass = (c1*b1 + c2*b2 + c3*j1 + c4*j2 + nu_corr + l).M();
        }
        else
        {
            tmp_hh_mass = (c1*b1 + c2*b2 + j1 + j2 + nu_corr + l).M();
        }
        hh_mass->Fill(tmp_hh_mass);
    }

    if (fails_counter != 0)
    {
        random_hme_failed = true;
        random_hme_eff = 1 - static_cast<float>(fails_counter)/nIter;
    }
    fails_counter = 0;

    int binmax = hh_mass->GetMaximumBin(); 
    float evt_hh_mass = hh_mass->GetXaxis()->GetBinCenter(binmax);

    return evt_hh_mass;
}

float hme::hme_simpl_impr(std::vector<TLorentzVector> const& particles, TH1F* h_mass, int nIter, TRandom3& rg)
{
    TLorentzVector b1 = particles[0];
    TLorentzVector b2 = particles[1];
    TLorentzVector j1 = particles[2];
    TLorentzVector j2 = particles[3];
    TLorentzVector l = particles[4];
    TLorentzVector met = particles[5];

    std::unique_ptr<TH1F> hh_mass = std::make_unique<TH1F>("evt_hh_mass", "evt_hh_mass", 80, 0.0, 2000.0);

    static int fails_counter = 0;

    rg.SetSeed(0);
    for (size_t i = 0; i < nIter; ++i)
    {
        float mh = h_mass->GetRandom();
        float universal_corr = mh/(b1+b2).M();
        b1 *= universal_corr;
        b2 *= universal_corr;

        float met_sigma = 25.2;

        float dpx = rg.Gaus(0, met_sigma);
        float dpy = rg.Gaus(0, met_sigma);

        float c1 = universal_corr;
        float c2 = universal_corr;

        float met_px_corr = -(c1 - 1)*b1.Px() - (c2 - 1)*b2.Px();
        float met_py_corr = -(c1 - 1)*b1.Py() - (c2 - 1)*b2.Py();

        met.SetPxPyPzE(met.Px() + dpx + met_px_corr, met.Py() + dpy + met_py_corr, met.Pz(), met.E());

        TLorentzVector vis = (l + j1 + j2);
        float a = mh*mh - vis.M()*vis.M() + 2.0*vis.Px()*met.Px() + 2.0*vis.Py()*met.Py();
        float A = 4.0*(vis.E()*vis.E() - vis.Pz()*vis.Pz());
        float B = -4.0*a*vis.Pz();
        float C = -4.0*vis.E()*vis.E()*(met.Px()*met.Px() + met.Py()*met.Py()) - a*a;
        float delta = B*B - 4.0*A*C;

        if (delta < 0.0f) 
        {
            ++fails_counter;
            continue;
        }

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

        float tmp_hh_mass = (l + nu + j1 + j2 + b1 + b2).M();
        hh_mass->Fill(tmp_hh_mass);
    }

    if (fails_counter != 0)
    {
        simpl_impr_hme_failed = true;
        simpl_impr_hme_eff = 1 - static_cast<float>(fails_counter)/nIter;
    }

    fails_counter = 0;

    int binmax = hh_mass->GetMaximumBin(); 
    float evt_hh_mass = hh_mass->GetXaxis()->GetBinCenter(binmax);

    return evt_hh_mass;
}