#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <optional>
#include <stdexcept>

#include "TRandom3.h"

#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"
#include "Constants.hpp"


class EstimatorBase
{
    public:
    EstimatorBase();

    private:
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass; 
};


class EstimatorSingleLep : public EstimatorBase
{
    public:
    EstimatorSingleLep(TString const& pdf_file_name);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets,
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb);

    // this method does not solve any constraints
    // it only assigns weights to each assignment of sampled parameters 
    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaWeights(VecLVF_t const& particles, 
                                                            std::pair<Float_t, Float_t> lj_pt_res, 
                                                            ULong64_t evt, 
                                                            TString const& comb_id);

    // this method computes an estimate of mass by solving physics consrtaints
    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaEqns(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
};


class EstimatorDoubleLep : public EstimatorBase
{
    public:
    EstimatorDoubleLep(TString const& pdf_file_name);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets,
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb);

    // this method does not solve any constraints
    // it only assigns weights to each assignment of sampled parameters 
    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaWeights(VecLVF_t const& particles, 
                                                            std::pair<Float_t, Float_t> lj_pt_res, 
                                                            ULong64_t evt, 
                                                            TString const& comb_id);

    // this method computes an estimate of mass by solving physics consrtaints
    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaEqns(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
};


class Estimator
{
    public:
    Estimator(TString const& pdf_file_name_sl, TString const& pdf_file_name_dl);

    private:
    EstimatorSingleLep m_estimator_sl;
    EstimatorDoubleLep m_estimator_dl;
};


class EstimatorSingLep_Run3
{
    public:
    EstimatorSingLep_Run3() = default;
    explicit EstimatorSingLep_Run3(TString const& file_name);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, std::vector<Float_t> const& jet_resolutions, LorentzVectorF_t const& met, ULong64_t evt, TString& chosen_comb);
    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, std::pair<Float_t, Float_t> lj_pt_res, ULong64_t evt, TString const& comb_id);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass;
};

class EstimatorDoubleLep_Run2
{
    public:
    EstimatorDoubleLep_Run2(TString const& pdf_file_name);

    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id);

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt,
                                        TString& chosen_comb);

    private:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass;
};

#endif