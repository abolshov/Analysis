#include "HistManager.hpp"

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

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
        auto&& [hist, hist_info] = m_hists.at(hist_name); 
        hist->Fill(value, weight);
    }
    catch(...)
    {
        return;
    }
}

void HistManager::Draw() const
{
    for (auto&& [hist_name, hhi_pair]: m_hists)
    {
        auto&& [hist, hist_info] = hhi_pair;
        auto&& [title, labels, range, nbins] = hist_info;
        auto&& [xlabel, ylabel] = labels;
        auto c1 = std::make_unique<TCanvas>("c1", "c1");
        c1->SetGrid();
        c1->SetTickx();
        c1->SetTicky();

        hist->GetXaxis()->SetTitle(xlabel.c_str());
        hist->GetYaxis()->SetTitle(ylabel.c_str());
        hist->SetLineWidth(3);
        hist->SetStats(1);
        hist->DrawNormalized();
        std::string fname = m_path;
        fname += "/";
        fname += hist_name;
        fname += ".png";
        c1->SaveAs(fname.c_str());
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
            auto&& [hist, hist_info] = m_hists.at(name); 
            auto&& [title, labels, range, nbins] = hist_info;
            auto&& [xlabel, ylabel] = labels;

            hist->SetLineWidth(3);
            hist->GetXaxis()->SetTitle(xlabel.c_str());
            hist->GetYaxis()->SetTitle(ylabel.c_str());
            if (color == 5) ++color;
            hist->SetLineColor(color++);
            stack->Add(hist.get());
            leg->AddEntry(hist.get(), name.c_str(), "l");
        }
        catch(...)
        {
            continue;
        }
    }

    stack->Draw("nostack");
    leg->Draw();
    c1->SaveAs((m_path + "/" + name + ".png").c_str());
}

void HistManager::DrawStack(std::vector<std::string>&& names, std::string&& title, std::string&& name) const
{
    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();  

    int color = 2;
    auto stack = std::make_unique<THStack>("stack", std::move(title).c_str());
    auto leg = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9); 

    for (auto&& name: names)
    { 
        try
        {
            auto&& [hist, hist_info] = m_hists.at(name); 
            auto&& [title, labels, range, nbins] = hist_info;
            auto&& [xlabel, ylabel] = labels;

            hist->SetLineWidth(3);
            hist->GetXaxis()->SetTitle(xlabel.c_str());
            hist->GetYaxis()->SetTitle(ylabel.c_str());
            if (color == 5) ++color;
            hist->SetLineColor(color++);
            stack->Add(hist.get());
            leg->AddEntry(hist.get(), name.c_str(), "l");
        }
        catch(...)
        {
            continue;
        }
    }

    stack->Draw("nostack");
    leg->Draw();
    c1->SaveAs((m_path + "/" + std::move(name) + ".png").c_str());
}