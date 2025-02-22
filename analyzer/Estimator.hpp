#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <optional>
#include <stdexcept>

#include "TRandom3.h"

#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"
#include "Constants.hpp"
#include "EstimationRecorder.hpp"


class EstimatorBase
{
    public:
    EstimatorBase();
    virtual ~EstimatorBase() = default;

    // this method does not solve any constraints
    // it only assigns weights to each assignment of sampled parameters 
    virtual std::array<Float_t, OUTPUT_SIZE> EstimateCombViaWeights(VecLVF_t const& particles, 
                                                                    std::pair<Float_t, Float_t> lj_pt_res, 
                                                                    ULong64_t evt, 
                                                                    TString const& comb_id) = 0;

    // this method computes an estimate of mass by solving physics consrtaints
    virtual std::array<Float_t, OUTPUT_SIZE> EstimateCombViaEqns(VecLVF_t const& particles, 
                                                                 ULong64_t evt, 
                                                                 TString const& comb_id) = 0;

    virtual std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                                VecLVF_t const& leptons, 
                                                std::vector<Float_t> const& jet_resolutions, 
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
    EstimationRecorder m_recorder; 

    virtual std::unique_ptr<TTree> MakeTree(TString const& tree_name) = 0;
};


class EstimatorSingleLep final : public EstimatorBase
{
    public:
    EstimatorSingleLep(TString const& pdf_file_name);

    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaWeights(VecLVF_t const& particles, 
                                                            std::pair<Float_t, Float_t> lj_pt_res, 
                                                            ULong64_t evt, 
                                                            TString const& comb_id) override;

    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaEqns(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id) override;

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> const& jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb) override;

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb) override;

    private: 
    std::array<LorentzVectorF_t, CONTROL_SL> j1{};
    std::array<LorentzVectorF_t, CONTROL_SL> j2{};
    std::array<LorentzVectorF_t, CONTROL_SL> met_corr{};
    std::array<LorentzVectorF_t, CONTROL_SL> nu{};
    std::array<LorentzVectorF_t, CONTROL_SL> lepW{};
    std::array<LorentzVectorF_t, CONTROL_SL> hadW{};
    std::array<LorentzVectorF_t, CONTROL_SL> Hww{};
    std::array<LorentzVectorF_t, CONTROL_SL> Xhh{};
    
    LorentzVectorF_t b1{};
    LorentzVectorF_t b2{};
    LorentzVectorF_t Hbb{};

    Float_t c3[CONTROL_SL] = {};
    Float_t c4[CONTROL_SL] = {};
    Float_t ljet_resc_dpx[CONTROL_SL] = {};
    Float_t ljet_resc_dpy[CONTROL_SL] = {};
    Float_t mass[CONTROL_SL] = {};

    Float_t c1{0.0};
    Float_t c2{0.0};
    Float_t mh{0.0};
    Float_t mw1{0.0};
    Float_t mw2{0.0};
    Float_t smear_dpx{0.0};
    Float_t smear_dpy{0.0};
    Float_t bjet_resc_dpx{0.0};
    Float_t bjet_resc_dpy{0.0};
    Float_t weight{0.0};
    Int_t num_sol{0};

    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;
};


class EstimatorDoubleLep final : public EstimatorBase
{
    public:
    EstimatorDoubleLep(TString const& pdf_file_name);

    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaWeights(VecLVF_t const& particles, 
                                                            std::pair<Float_t, Float_t> lj_pt_res, 
                                                            ULong64_t evt, 
                                                            TString const& comb_id) override;

    std::array<Float_t, OUTPUT_SIZE> EstimateCombViaEqns(VecLVF_t const& particles, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id) override;

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> const& jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb) override;

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb) override;

    private: 
    std::array<LorentzVectorF_t, CONTROL_DL> nu_offshell{};
    std::array<LorentzVectorF_t, CONTROL_DL> nu_onshell{};
    std::array<LorentzVectorF_t, CONTROL_DL> l_offshell{};
    std::array<LorentzVectorF_t, CONTROL_DL> l_onshell{};
    std::array<LorentzVectorF_t, CONTROL_DL> offshellW{};
    std::array<LorentzVectorF_t, CONTROL_DL> onshellW{};
    std::array<LorentzVectorF_t, CONTROL_DL> Hww{};
    std::array<LorentzVectorF_t, CONTROL_DL> Hbb{};
    std::array<LorentzVectorF_t, CONTROL_DL> Xhh{};
    
    LorentzVectorF_t b1{};
    LorentzVectorF_t b2{};
    LorentzVectorF_t met_corr{};

    Float_t c1{0.0};
    Float_t c2{0.0};
    Float_t mw{0.0};
    Float_t mh{0.0};
    Float_t phi_gen{0.0};
    Float_t eta_gen{0.0};
    Float_t smear_dpx{0.0};
    Float_t smear_dpy{0.0};
    Float_t bjet_resc_dpx{0.0};
    Float_t bjet_resc_dpy{0.0};
    Float_t weight{0.0};
    Float_t mass{0.0};
    Int_t num_sol{0};

    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;
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

    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, 
                                        VecLVF_t const& leptons, 
                                        std::vector<Float_t> const& jet_resolutions, 
                                        LorentzVectorF_t const& met, 
                                        ULong64_t evt, 
                                        TString& chosen_comb);

    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, 
                                                         std::pair<Float_t, Float_t> lj_pt_res, 
                                                         ULong64_t evt, 
                                                         TString const& comb_id);

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