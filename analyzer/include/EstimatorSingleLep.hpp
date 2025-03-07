#ifndef ESTIM_SL
#define ESTIM_SL

#include "EstimatorBase.hpp"
#include "Constants.hpp"
#include "Definitions.hpp"

class EstimatorSingleLep final : public EstimatorBase
{
    public:
    EstimatorSingleLep(TString const& pdf_file_name, TString const& dbg_file_name);
    ~EstimatorSingleLep() override = default;

    std::array<Float_t, OUTPUT_SIZE> EstimateCombination(VecLVF_t const& particles, ULong64_t evt_id, TString const& comb_label) override;
    std::optional<Float_t> EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id) override;

    private: 
    struct IterData;
    std::unique_ptr<IterData> m_iter_data;
    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;
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

#endif