#include "HistManager.hpp"

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"

HistManager::HistManager(int n) 
{
    m_hists.reserve(n); 
}

void HistManager::Add(std::string const& hist_name, std::string const& title, std::pair<std::string, std::string> labels, std::pair<double, double> range, int n_bins)
{
    HistInfo hi{title, labels, range, n_bins};
    auto [min, max] = range;
    m_hists.emplace(hist_name, std::make_pair(std::make_unique<TH1F>(hist_name.data(), title.data(), n_bins, min, max), hi));
}

void HistManager::Fill(std::string const& hist_name, double value, double weight) noexcept
{
    try
    {
        auto const& [hist, hist_info] = m_hists.at(hist_name); 
        hist->Fill(value, weight);
    }
    catch(...)
    {
        return;
    }
}

void HistManager::Draw() const
{
    for (auto const& [hist_name, hhi_pair]: m_hists)
    {
        auto const& [ptr_hist, hist_info] = hhi_pair;
        auto const& [title, labels, range, nbins] = hist_info;
        auto const& [xlabel, ylabel] = labels;
        auto c1 = std::make_unique<TCanvas>("c1", "c1");
        c1->SetGrid();
        c1->SetTickx();
        c1->SetTicky();

        ptr_hist->GetXaxis()->SetTitle(xlabel.c_str());
        ptr_hist->GetYaxis()->SetTitle(ylabel.c_str());
        ptr_hist->SetLineWidth(3);
        ptr_hist->SetStats(1);
        ptr_hist->DrawNormalized();
        c1->SaveAs(Form("histograms/%s.png", hist_name.c_str()));
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

    for (auto const& name: names)
    { 
        try
        {
            auto const& [ptr_hist, hist_info] = m_hists.at(name); 
            auto const& [title, labels, range, nbins] = hist_info;
            auto const& [xlabel, ylabel] = labels;

            ptr_hist->SetLineWidth(3);
            ptr_hist->GetXaxis()->SetTitle(xlabel.c_str());
            ptr_hist->GetYaxis()->SetTitle(ylabel.c_str());
            if (color == 5) ++color;
            ptr_hist->SetLineColor(color++);
            stack->Add(ptr_hist.get());
            leg->AddEntry(ptr_hist.get(), name.c_str(), "l");
        }
        catch(...)
        {
            continue;
        }
    }

    stack->Draw("nostack");
    leg->Draw();
    c1->SaveAs(Form("histograms/%s.png", name.c_str()));
}