#include "HistManager.hpp"

#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TFile.h"

void HistManager::Add(std::string hist_name, std::string const& title, Label const& labels, Range const& range, int n_bins)
{
    m_hists_1d.insert({std::move(hist_name), MakeItem1D(hist_name, title, labels, range, n_bins)});
}

void HistManager::Add(std::string hist_name, std::string const& title, Label const& labels, Range const& xrange, Range const& yrange, Bins const& bins)
{
    m_hists_2d.insert({std::move(hist_name), MakeItem2D(hist_name, title, labels, xrange, yrange, bins)});
}

void HistManager::Fill(std::string const& hist_name, double value)
{
    auto const& [hist, hist_info] = m_hists_1d.at(hist_name); 
    hist->Fill(value);
}

void HistManager::Fill(std::string const& hist_name, double xval, double yval)
{
    auto const& [hist, hist_info] = m_hists_2d.at(hist_name); 
    hist->Fill(xval, yval);
}

void HistManager::FillWeighted(std::string const& hist_name, double value, double weight)
{
    auto const& [hist, hist_info] = m_hists_1d.at(hist_name); 
    hist->Fill(value, weight);
}

void HistManager::FillWeighted(std::string const& hist_name, double xval, double yval, double weight)
{
    auto const& [hist, hist_info] = m_hists_2d.at(hist_name); 
    hist->Fill(xval, yval, weight);
}

void HistManager::Draw() const
{
    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);

    for (auto const& [hist_name, hhi_pair]: m_hists_1d)
    {
        auto const& [ptr_hist, hist_info] = hhi_pair;
        auto const& [title, labels, range, nbins] = hist_info;
        auto const& [xlabel, ylabel] = labels;

        ptr_hist->GetXaxis()->SetTitle(xlabel.c_str());
        ptr_hist->GetYaxis()->SetTitle(ylabel.c_str());
        ptr_hist->SetLineWidth(3);
        ptr_hist->SetStats(1);
        ptr_hist->Draw("hist");
        // ptr_hist->DrawNormalized();
        c1->SaveAs(Form("histograms/%s.png", hist_name.c_str()));
        // c1->SaveAs(Form("%s/%s.png", m_path.c_str(), hist_name.c_str()));
    }

    for (auto const& [hist_name, hhi_pair]: m_hists_2d)
    {
        auto const& [ptr_hist, hist_info] = hhi_pair;
        auto const& [title, labels, xrange, yrange, bins] = hist_info;
        auto const& [xlabel, ylabel] = labels;
        // auto c1 = std::make_unique<TCanvas>("c1", "c1");
        // c1->SetGrid();
        // c1->SetTickx();
        // c1->SetTicky();

        ptr_hist->GetXaxis()->SetTitle(xlabel.c_str());
        ptr_hist->GetYaxis()->SetTitle(ylabel.c_str());
        ptr_hist->SetStats(1);
        ptr_hist->Draw("colz");
        c1->SaveAs(Form("histograms/%s.png", hist_name.c_str()));
        // c1->SaveAs(Form("%s/%s.png", m_path.c_str(), hist_name.c_str()));
    }
}

void HistManager::DrawStack(std::vector<std::string> const& names, std::string const& title, std::string const& name) const
{
    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();  

    int color = 2;
    auto stack = std::make_unique<THStack>("stack", title.c_str());
    auto leg = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);

    std::string stack_xlabel, stack_ylabel;

    for (auto const& name: names)
    { 
        try
        {
            auto const& [ptr_hist, hist_info] = m_hists_1d.at(name); 
            auto const& [title, labels, range, nbins] = hist_info;
            auto const& [xlabel, ylabel] = labels;

            stack_xlabel =  xlabel; 
            stack_ylabel =  ylabel; 

            ptr_hist->SetLineWidth(3);
            ptr_hist->SetStats(0);
            if (color == 5) ++color;
            ptr_hist->SetLineColor(color++);
            stack->Add(ptr_hist.get(), "hist");
            leg->AddEntry(ptr_hist.get(), name.c_str(), "l");
        }
        catch(...)
        {
            continue;
        }
    }
    
    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle(stack_xlabel.c_str());
    stack->GetYaxis()->SetTitle(stack_ylabel.c_str());
    leg->Draw();
    c1->SaveAs(Form("histograms/%s.png", name.c_str()));
    // c1->SaveAs(Form("%s/%s.png", m_path.c_str(), name.c_str()));
}

void HistManager::Write1D(std::string const& fname, std::vector<std::string> const& names) const
{
    auto output = std::make_unique<TFile>(fname.c_str(), "RECREATE");
    for (auto const& name: names)
    {
        auto const& [ptr_hist, hist_info] = m_hists_1d.at(name);
        int binmax = ptr_hist->GetMaximumBin();
        double content = ptr_hist->GetBinContent(binmax);
        ptr_hist->Scale(1.0/content);
        ptr_hist->Write();
    }
	output->Write();
	output->Close();
}

void HistManager::Write2D(std::string const& fname, std::vector<std::string> const& names) const
{
    auto output = std::make_unique<TFile>(fname.c_str(), "UPDATE");
    for (auto const& name: names)
    {
        auto const& [ptr_hist, hist_info] = m_hists_2d.at(name);
        ptr_hist->Write();
    }
	output->Write();
	output->Close();
}