#ifndef ESTIMATOR_UTILS_HPP
#define ESTIMATOR_UTILS_HPP

#include <vector>
#include <type_traits>

#include "Definitions.hpp"
#include "Storage.hpp"
#include "Constants.hpp"

#include "TFile.h"
#include "TH1.h"

template <typename T>
UHist_t<T> ReadHist(TFile* file, TString hist_name)
{
    UHist_t<T> res(static_cast<T*>(file->Get<T>(hist_name)));
    res->SetDirectory(nullptr);
    return res;
}

template <typename T>
UHist_t<T> Copy(UHist_t<T> const& hist)
{
    UHist_t<T> res(static_cast<T*>(hist->Clone()));
    res->SetDirectory(nullptr);
    return res;
}

template <typename T>
void ResetHist(UHist_t<T>& hist)
{
    hist->Reset("ICES");
}

void Get1dPDFs(TFile* fptr, HistVec_t<TH1F>& pdfs, Channel ch);
void Get2dPDFs(TFile* fptr, HistVec_t<TH2F>& pdfs, Channel ch);

inline LorentzVectorF_t GetRecoMET(Storage const& s) { return LorentzVectorF_t(s.reco_met_pt, 0.0, s.reco_met_phi, 0.0); }
VecLVF_t GetRecoJetP4(Storage const& s);
VecLVF_t GetRecoLepP4(Storage const& s, Channel ch);
std::vector<Float_t> GetPNetRes(Storage const& s);

#ifdef DEV
    // inline LorentzVectorF_t GetGenMET(Storage const& s) { return LorentzVectorF_t(s.gen_met_pt, 0.0, s.gen_met_phi, 0.0); }
    // VecLVF_t GetGenJetP4(Storage const& s);
    VecLVF_t GetGenLepP4(Storage const& s, Channel ch);
    VecLVF_t GetGenQuarksP4(Storage const& s, Channel ch);
    VecLVF_t GetGenNuP4(Storage const& s, Channel ch);
#endif

Float_t ComputeWidth(UHist_t<TH1F> const& h, unsigned l, unsigned r);

template <typename It, std::enable_if_t<std::is_floating_point_v<typename It::value_type>, bool> = true>
void ZScoreTransform(It begin, It end)
{
    if (begin == end)
    {
        return;
    }

    typename It::value_type sum = 0.0;
    typename It::value_type sum_sqr = 0.0;

    It it = begin;
    while (it != end)
    {
        typename It::value_type val = *it;
        sum += val;
        sum_sqr += val*val;
        ++it;
    }

    size_t n = end - begin;
    typename It::value_type mean = sum/n;
    typename It::value_type mean_sqr = sum_sqr/n;
    typename It::value_type std = std::sqrt(mean_sqr - mean*mean);

    it = begin;
    while (it != end)
    {
        *it -= mean;
        *it /= std;
        ++it;
    }
}

template <typename It, std::enable_if_t<std::is_floating_point_v<typename It::value_type>, bool> = true>
void MinMaxTransform(It begin, It end)
{
    if (begin == end)
    {
        return;
    }

    auto [min_it, max_it] = std::minmax_element(begin, end);
    typename It::value_type diff = *max_it - *min_it;
    typename It::value_type min = *min_it;

    auto Func = [&min, &diff](typename It::value_type const& val)
    { 
        typename It::value_type ret = val;
        ret -= min; 
        ret /= diff;
        return ret;
    };
    std::transform(begin, end, begin, Func);
}
#endif