#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <optional>
#include <stdexcept>

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

    virtual std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                                VecLVF_t const& leptons, 
                                                std::vector<Float_t> jet_resolutions, 
                                                LorentzVectorF_t const& met, 
                                                ULong64_t evt,
                                                TString& chosen_comb) = 0;

    virtual std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                                VecLVF_t const& leptons, 
                                                LorentzVectorF_t const& met, 
                                                ULong64_t evt,
                                                TString& chosen_comb) = 0;

    protected:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass; 
};

// class EstimatorSingLep_Run3 final : public Estimator
class EstimatorSingLep_Run3
{
    public:
    EstimatorSingLep_Run3() = default;
    explicit EstimatorSingLep_Run3(TString const& file_name);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, std::vector<Float_t> jet_resolutions, LorentzVectorF_t const& met, ULong64_t evt, TString& chosen_comb);
    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, std::pair<Float_t, Float_t> lj_pt_res, ULong64_t evt, TString const& comb_id);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass;

    // std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt) override
    // {
    //     std::cout << "Illegal to call this function here\n";
    // }
};

class EstimatorDoubleLep_Run2 final : public Estimator
{
    public:
    // EstimatorDoubleLep_Run2() = default;
    using Estimator::Estimator;

    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt,
                                        TString& chosen_comb) override;

    private:
    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt,
                                        TString& chosen_comb) override
    {
        throw std::runtime_error("Illegal function call");
    }
};

#endif