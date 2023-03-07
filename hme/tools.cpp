#include "tools.hpp"
#include <iostream>

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
    
    float mass;

    for (size_t i = 0; i < 1000; ++i)
    {
        float lead_rescale_factor = lead_jet_pdf->GetRandom();
        
        mass = mass_pdf->GetRandom(); 

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

    mass = mass_pdf->GetRandom();

    float universal_corr = mass/(p1+p2).M();
    return std::make_pair(universal_corr, universal_corr);

    // if (x2 < 0)
    // {
    //     std::cout << "No solution: negative scalar product of jet momenta" << std::endl; 
    //     return std::make_pair(universal_corr, universal_corr);
    // }
    // if ((x2*x2 - 4*x1*x3) < 0 || x1 == 0)
    // {
    //     std::cout << "No solution: negative discriminiant" << std::endl;
    //     return std::make_pair(universal_corr, universal_corr);
    // }
    // float trail_rescale_factor = (-x2 + sqrt(x2*x2 - 4*x1*x3))/(2*x1);
    // if (trail_rescale_factor < 0)
    // {
    //     std::cout << "No solution: negative trailing jet rescaling factor" << std::endl;
    //     return std::make_pair(universal_corr, universal_corr);
    // }
    // return std::make_pair(lead_rescale_factor, trail_rescale_factor);
}