#ifndef COMB_TOOLS_HPP
#define COMB_TOOLS_HPP

#include <vector>
#include <algorithm>
#include <numeric>

template <typename T>
std::vector<int> sort_indices(T* v, int n) 
{
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2];});
  return idx;
}

template <typename Func>
std::pair<int, int> ChooseBestPair(std::vector<TLorentzVector> const& jets, Func func)
{
    size_t sz = jets.size();
    if (sz == 2)
    {
        return {0, 1};
    }


    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(sz*(sz-1)/2);
    for (size_t i = 0; i < sz; ++i)
    {
        for (size_t j = 0; j < sz; ++j)
        {
            pairs.emplace_back(i, j);
        }
    }

    double min_metric = 10e9;
    std::pair<int, int> res{-1, -1};
    for (auto const& p: pairs)
    {
        auto [i1, i2] = p;
        TLorentzVector const& jet1 = jets[i1];
        TLorentzVector const& jet2 = jets[i2];

        double metric = func(jet1, jet2);
        if (metric < min_metric)
        {
            min_metric = metric;
            res = p;
        }
    }

    return res;
}

std::pair<int, int> FindByAngle(std::vector<TLorentzVector> const& jets, int bj1_idx, int bj2_idx)
{
    TLorentzVector Hbb = jets[bj1_idx] + jets[bj2_idx];
    int lj1_idx = -1;
    int lj2_idx = -1;

    double max_dphi = -5.0;

    int sz = jets.size();
    for (int i = 0; i < sz; ++i)
    {
        if (i == bj1_idx || i == bj2_idx)
        {
            continue;
        }

        for (int j = i + 1; j < sz; ++j)
        {
            if (j == bj1_idx || j == bj2_idx)
            {
                continue;
            }

            TLorentzVector Hww = jets[i] + jets[j];
            double dphi = std::abs(Hww.DeltaPhi(Hbb));
            if (dphi > max_dphi)
            {
                max_dphi = dphi;
                lj1_idx = i;
                lj2_idx = j;
            }
        }
    }

    return {lj1_idx, lj2_idx};
}

// make all possible combinations of light jet pairs given indices of b jets
std::vector<std::pair<int, int>> MakePairs(int sz, int skip1, int skip2)
{
    std::vector<std::pair<int, int>> pairs;
    for (int i = 0; i < sz; ++i)
    {
        if (i == skip1 || i == skip2)
        {
            continue;
        }

        for (int j = i + 1; j < sz; ++j)
        {
            if (j == skip1 || j == skip2)
            {
                continue;
            }
            pairs.emplace_back(i, j);
        }
    }
    return pairs;
}

#endif