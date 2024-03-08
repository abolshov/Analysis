#ifndef HIST_MAN_HPP
#define HIST_MAN_HPP

#include <unordered_map>
#include <memory>

#include "TH1.h"
#include "TH2.h"

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

    using Item1D = std::pair<std::unique_ptr<TH1F>, HistInfo>;
    using Item2D = std::pair<std::unique_ptr<TH2F>, HistInfo>;

    inline HistInfo MakeInfo(std::string const& title, std::pair<std::string, std::string> const& labels, std::pair<double, double> const& range, int n_bins)
    {
        return {title, labels, range, n_bins};
    }

    inline std::unique_ptr<TH1F> Make1DHist(std::string const& name, std::string const& title, std::pair<double, double> const& range, int n_bins)
    {
        auto [min, max] = range;
        return std::make_unique<TH1F>(name.data(), title.data(), n_bins, min, max);
    }

    inline Item1D MakeItem1D(std::string const& hist_name, std::string const& title, std::pair<std::string, std::string> const& labels, std::pair<double, double> const& range, int n_bins)
    {
        return {Make1DHist(hist_name, title, range, n_bins), MakeInfo(title, labels, range, n_bins)};
    }

    std::unordered_map<std::string, Item1D> m_hists_1d;
    // std::unordered_map<std::string, Item2D> m_hists_2d;

    public:
    void Add(std::string hist_name, std::string const& title, std::pair<std::string, std::string> const& labels, std::pair<double, double> const& range, int n_bins);
    void Fill(std::string const& hist_name, double value, double weight = 1.0);
    void Draw() const;
    void DrawStack(std::vector<std::string> const& names, std::string const& title, std::string const& name) const;
};

#endif