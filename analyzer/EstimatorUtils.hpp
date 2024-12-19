#ifndef ESTIMATOR_UTILS_HPP
#define ESTIMATOR_UTILS_HPP

#include <vector>
#include <type_traits>

#include "Definitions.hpp"
#include "Storage.hpp"

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
void Reset(UHist_t<T>& hist)
{
    hist->Reset("ICES");
}

void Get1dPDFs(TFile* fptr, HistVec_t<TH1F>& pdfs);
void Get2dPDFs(TFile* fptr, HistVec_t<TH2F>& pdfs);

inline LorentzVectorF_t GetRecoMET(Storage const& s) { return LorentzVectorF_t(s.reco_met_pt, 0.0, s.reco_met_phi, 0.0); }
VecLVF_t GetRecoJetP4(Storage const& s);
VecLVF_t GetRecoLepP4(Storage const& s, Channel ch);
std::vector<Float_t> GetPNetRes(Storage const& s);

#ifdef DEBUG
inline LorentzVectorF_t GetGenMET(Storage const& s) { return LorentzVectorF_t(s.gen_met_pt, 0.0, s.gen_met_phi, 0.0); }
VecLVF_t GetGenJetP4(Storage const& s);
VecLVF_t GetGenLepP4(Storage const& s, Channel ch);
VecLVF_t GetGenQuarksP4(Storage const& s, Channel ch);
VecLVF_t GetGenNuP4(Storage const& s, Channel ch);
#endif

Float_t ComputeWidth(UHist_t<TH1F> const& h, unsigned l, unsigned r);

template <typename It, std::enable_if_t<std::is_arithmetic_v<typename It::value_type>, bool> = true>
void ZScoreTransform(It begin, It end)
{
    if (begin == end)
    {
        return;
    }

    Float_t sum = 0.0;
    Float_t sum_sqr = 0.0;

    It it = begin;
    while (it != end)
    {
        Float_t val = static_cast<Float_t>(*it);
        sum += val;
        sum_sqr += val*val;
        ++it;
    }

    size_t n = end - begin;
    Float_t mean = sum/n;
    Float_t mean_sqr = sum_sqr/n;
    Float_t std = std::sqrt(mean_sqr - mean*mean);

    it = begin;
    while (it != end)
    {
        *it -= mean;
        *it /= std;
        ++it;
    }
}

template <typename It, std::enable_if_t<std::is_arithmetic_v<typename It::value_type>, bool> = true>
void MinMaxTransform(It begin, It end)
{
    if (begin == end)
    {
        return;
    }

    auto [min_it, max_it] = std::minmax_element(begin, end);
    auto diff = *max_it - *min_it;

    It it = begin;
    while (it != end)
    {
        *it = (*it - *min_it)/diff;
        ++it;
    }
}

#endif