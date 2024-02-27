#ifndef HIST_MAN_HPP
#define HIST_MAN_HPP

#include <unordered_map>
#include <memory>

#include "TH1.h"

class HistManager
{
    private:
    struct HistInfo
    {
        std::string m_title;
        std::pair<std::string, std::string> m_axis_labels;
        std::pair<double, double> m_range;
        int n_bins;
    };

    std::unordered_map<std::string, std::pair<std::unique_ptr<TH1F>, HistInfo>> m_hists;

    public:
    explicit HistManager(int n); 

    HistManager(HistManager const& other) = default;
    HistManager(HistManager&& other) = default;
    HistManager& operator=(HistManager const& other) = default;
    HistManager& operator=(HistManager&& other) = default;

    void Add(std::string const& hist_name, std::string const& title, std::pair<std::string, std::string> const& labels, std::pair<double, double> range, int n_bins);
    void Fill(std::string const& hist_name, double value, double weight = 1.0) noexcept;
    void Draw() const;
    void DrawStack(std::vector<std::string> const& names, std::string const& title, std::string const& name) const;
};

#endif