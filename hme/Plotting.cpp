#include "Plotting.hpp"

void save_1d_dist(TH1F* dist, std::string const& name, std::string const& title)
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

void save_2d_dist(TH2F* dist, std::string const& name, std::string const& title_x, std::string const& title_y)
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