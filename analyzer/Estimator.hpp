#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include <optional>
#include <type_traits>

#include "TRandom3.h"

#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "EstimatorTools.hpp"
#include "Constants.hpp"
#include "EstimationRecorder.hpp"


template <typename T, std::enable_if_t<std::is_default_constructible_v<T>, bool> = true>
void ResetObject(T& object)
{
    object = T{};
}


class EstimatorBase
{
    public:
    explicit EstimatorBase(TString dbg_file_name = {});
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
    bool m_record;

    virtual std::unique_ptr<TTree> MakeTree(TString const& tree_name) = 0;
};


class EstimatorSingleLep final : public EstimatorBase
{
    public:
    EstimatorSingleLep(TString const& pdf_file_name, TString const& dbg_file_name);

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
    struct IterData;
    std::unique_ptr<IterData> m_iter_data;

    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;
};


class EstimatorDoubleLep final : public EstimatorBase
{
    public:
    EstimatorDoubleLep(TString const& pdf_file_name, TString const& dbg_file_name);

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
    struct IterData;
    std::unique_ptr<IterData> m_iter_data;

    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;
};


class Estimator
{
    public:
    Estimator(TString const& pdf_file_name_sl, TString const& pdf_file_name_dl, TString const& dbg_file_name);

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

struct EstimatorSingleLep::IterData
{
    Float_t j1_pt[CONTROL] = {};
    Float_t j1_eta[CONTROL] = {};
    Float_t j1_phi[CONTROL] = {};
    Float_t j1_mass[CONTROL] = {};

    Float_t j2_pt[CONTROL] = {};
    Float_t j2_eta[CONTROL] = {};
    Float_t j2_phi[CONTROL] = {};
    Float_t j2_mass[CONTROL] = {};

    Float_t lepW_pt[CONTROL] = {};
    Float_t lepW_eta[CONTROL] = {};
    Float_t lepW_phi[CONTROL] = {};
    Float_t lepW_mass[CONTROL] = {};

    Float_t hadW_pt[CONTROL] = {};
    Float_t hadW_eta[CONTROL] = {};
    Float_t hadW_phi[CONTROL] = {};
    Float_t hadW_mass[CONTROL] = {};

    Float_t Hww_pt[CONTROL] = {};
    Float_t Hww_eta[CONTROL] = {};
    Float_t Hww_phi[CONTROL] = {};
    Float_t Hww_mass[CONTROL] = {};

    Float_t Xhh_pt[CONTROL] = {};
    Float_t Xhh_eta[CONTROL] = {};
    Float_t Xhh_phi[CONTROL] = {};
    Float_t Xhh_mass[CONTROL] = {};

    Float_t nu_pt[CONTROL] = {};
    Float_t nu_eta[CONTROL] = {};
    Float_t nu_phi[CONTROL] = {};

    Float_t mass[CONTROL] = {};

    Float_t met_corr_pt[CONTROL] = {};
    Float_t met_corr_phi[CONTROL] = {};

    Float_t ljet_resc_fact_1[CONTROL] = {};
    Float_t ljet_resc_fact_2[CONTROL] = {};

    Float_t ljet_resc_dpx[CONTROL] = {};
    Float_t ljet_resc_dpy[CONTROL] = {};

    Float_t b1_pt{0.0};
    Float_t b1_eta{0.0};
    Float_t b1_phi{0.0};
    Float_t b1_mass{0.0};

    Float_t b2_pt{0.0};
    Float_t b2_eta{0.0};
    Float_t b2_phi{0.0};
    Float_t b2_mass{0.0};

    Float_t Hbb_pt{0.0};
    Float_t Hbb_eta{0.0};
    Float_t Hbb_phi{0.0};
    Float_t Hbb_mass{0.0};

    Float_t bjet_resc_fact_1{0.0};
    Float_t bjet_resc_fact_2{0.0};
    Float_t mh{0.0};
    Float_t mw1{0.0};
    Float_t mw2{0.0};
    Float_t smear_dpx{0.0};
    Float_t smear_dpy{0.0};
    Float_t bjet_resc_dpx{0.0};
    Float_t bjet_resc_dpy{0.0};
    Float_t weight{0.0};
    Int_t num_sol{0};
};

struct EstimatorDoubleLep::IterData
{
    Float_t offshellW_pt[CONTROL] = {};
    Float_t offshellW_eta[CONTROL] = {};
    Float_t offshellW_phi[CONTROL] = {};
    Float_t offshellW_mass[CONTROL] = {};

    Float_t onshellW_pt[CONTROL] = {};
    Float_t onshellW_eta[CONTROL] = {};
    Float_t onshellW_phi[CONTROL] = {};
    Float_t onshellW_mass[CONTROL] = {};

    Float_t Hww_pt[CONTROL] = {};
    Float_t Hww_eta[CONTROL] = {};
    Float_t Hww_phi[CONTROL] = {};
    Float_t Hww_mass[CONTROL] = {};

    Float_t Xhh_pt[CONTROL] = {};
    Float_t Xhh_eta[CONTROL] = {};
    Float_t Xhh_phi[CONTROL] = {};
    Float_t Xhh_mass[CONTROL] = {};

    Float_t nu_offshell_pt[CONTROL] = {};
    Float_t nu_offshell_eta[CONTROL] = {};
    Float_t nu_offshell_phi[CONTROL] = {};

    Float_t nu_onshell_pt[CONTROL] = {};
    Float_t nu_onshell_eta[CONTROL] = {};
    Float_t nu_onshell_phi[CONTROL] = {};

    Float_t mass[CONTROL] = {};

    Float_t met_corr_pt{0.0};
    Float_t met_corr_phi{0.0};

    Float_t b1_pt{0.0};
    Float_t b1_eta{0.0};
    Float_t b1_phi{0.0};
    Float_t b1_mass{0.0};

    Float_t b2_pt{0.0};
    Float_t b2_eta{0.0};
    Float_t b2_phi{0.0};
    Float_t b2_mass{0.0};

    Float_t Hbb_pt{0.0};
    Float_t Hbb_eta{0.0};
    Float_t Hbb_phi{0.0};
    Float_t Hbb_mass{0.0};

    Float_t bjet_resc_fact_1{0.0};
    Float_t bjet_resc_fact_2{0.0};
    Float_t mh{0.0};
    Float_t mw{0.0};
    Float_t smear_dpx{0.0};
    Float_t smear_dpy{0.0};
    Float_t bjet_resc_dpx{0.0};
    Float_t bjet_resc_dpy{0.0};
    Float_t weight{0.0};
    Int_t num_sol{0};
};

#endif