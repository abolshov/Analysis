#ifndef HIST_MAN_HPP
#define HIST_MAN_HPP

#include <unordered_map>
#include <memory>
#include <type_traits>

#include "TH1.h"
#include "TH2.h"

class HistManager
{
    private:
    
    using Range = std::pair<double, double>;
    using Bins = std::pair<int, int>;
    using Label = std::pair<std::string, std::string>;

    struct Hist1DInfo
    {
        std::string m_title;
        Label m_axis_labels;
        Range m_range;
        int n_bins;
    };

    struct Hist2DInfo
    {
        std::string m_title;
        Label m_axis_labels;
        Range m_xrange;
        Range m_yrange;
        Bins m_bins;
    };

    using Item1D = std::pair<std::unique_ptr<TH1F>, Hist1DInfo>;
    using Item2D = std::pair<std::unique_ptr<TH2F>, Hist2DInfo>;

    inline Hist1DInfo Make1DInfo(std::string const& title, Label const& labels, Range const& range, int n_bins)
    {
        return {title, labels, range, n_bins};
    }

    inline Hist2DInfo Make2DInfo(std::string const& title, Label const& labels, Range const& xrange, Range const& yrange, Bins const& bins)
    {
        return {title, labels, xrange, yrange, bins};
    }

    inline std::unique_ptr<TH1F> Make1DHist(std::string const& name, std::string const& title, Range const& range, int n_bins)
    {
        auto const& [min, max] = range;
        return std::make_unique<TH1F>(name.data(), title.data(), n_bins, min, max);
    }

    inline std::unique_ptr<TH2F> Make2DHist(std::string const& name, std::string const& title, Range const& xrange, Range const& yrange, Bins const& bins)
    {
        auto const& [xmin, xmax] = xrange;
        auto const& [ymin, ymax] = yrange;
        auto const& [xbins, ybins] = bins;
        return std::make_unique<TH2F>(name.data(), title.data(), xbins, xmin, xmax, ybins, ymin, ymax);
    }

    inline Item1D MakeItem1D(std::string const& hist_name, std::string const& title, Label const& labels, Range const& range, int n_bins)
    {
        return {Make1DHist(hist_name, title, range, n_bins), Make1DInfo(title, labels, range, n_bins)};
    }

    inline Item2D MakeItem2D(std::string const& hist_name, std::string const& title, Label const& labels, Range const& xrange, Range const& yrange, Bins const& bins)
    {
        return {Make2DHist(hist_name, title, xrange, yrange, bins), Make2DInfo(title, labels, xrange, yrange, bins)};
    }

    std::unordered_map<std::string, Item1D> m_hists_1d;
    std::unordered_map<std::string, Item2D> m_hists_2d;

    public:
    void Add(std::string hist_name, std::string const& title, Label const& labels, Range const& range, int n_bins);
    void Add(std::string hist_name, std::string const& title, Label const& labels, Range const& xrange, Range const& yrange, Bins const& bins);
    void Fill(std::string const& hist_name, double value);
    void FillWeighted(std::string const& hist_name, double value, double weight);
    void Fill(std::string const& hist_name, double xval, double yval);
    void FillWeighted(std::string const& hist_name, double xval, double yval, double weight);
    void Draw() const;
    void DrawStack(std::vector<std::string> const& names, std::string const& title, std::string const& name) const;
    
    template <typename T>
    std::unique_ptr<T> const& Get(std::string const& name) const
    {
        if constexpr (std::is_same_v<T, TH1F>)
        {
            return m_hists_1d.at(name).first;
        }
        if constexpr (std::is_same_v<T, TH2F>)
        {
            return m_hists_2d.at(name).first;
        }
    }
};

#endif