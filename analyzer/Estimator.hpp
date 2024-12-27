#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <optional>

#include "TRandom3.h"

#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"
#include "Constants.hpp"


class Estimator
{
    public:
    Estimator() = default;
    explicit Estimator(TString const& file_name);
    virtual ~Estimator() = default;

    virtual std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, std::vector<Float_t> jet_resolutions, LorentzVectorF_t const& met, ULong64_t evt) = 0;

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass; 
};

class EstimatorSingLep
{
    public:
    EstimatorSingLep() = default;
    explicit EstimatorSingLep(TString const& file_name);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, std::vector<Float_t> jet_resolutions, LorentzVectorF_t const& met, ULong64_t evt, TString& chosen_comb);
    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, std::pair<Float_t, Float_t> lj_pt_res, ULong64_t evt, TString const& comb_id);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass;
};

#endif