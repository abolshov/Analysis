#ifndef ESTIM_DL
#define ESTIM_DL

#include "EstimatorBase.hpp"
#include "Constants.hpp"
#include "Definitions.hpp"

class EstimatorDoubleLep final : public EstimatorBase
{
    public:
    EstimatorDoubleLep(TString const& pdf_file_name, AggregationMode aggr_mode);
    ~EstimatorDoubleLep() override = default;

    OptArrF_t<ESTIM_OUT_SZ> EstimateMass(Event const& event) override;

    private: 
    // fix: write variables into this struct in EstimateCombSlim and EstimateCombFat directly instead of local vars
    struct IterData;
    std::unique_ptr<IterData> m_iter_data; 
    std::unique_ptr<TTree> MakeTree(TString const& tree_name) override;

    ArrF_t<ESTIM_OUT_SZ> EstimateCombSlim(VecLVF_t const& particles, 
                                          [[maybe_unused]] std::vector<Float_t> const& jet_res, 
                                          ULong64_t evt_id, 
                                          JetComb const& comb);

    ArrF_t<ESTIM_OUT_SZ> EstimateCombFat(VecLVF_t const& particles, 
                                         [[maybe_unused]] std::vector<Float_t> const& jet_res, 
                                         ULong64_t evt_id, 
                                         [[maybe_unused]] JetComb const& comb);
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

    Float_t l_offshell_pt[CONTROL] = {};
    Float_t l_offshell_eta[CONTROL] = {};
    Float_t l_offshell_phi[CONTROL] = {};
    Float_t l_offshell_mass[CONTROL] = {};

    Float_t l_onshell_pt[CONTROL] = {};
    Float_t l_onshell_eta[CONTROL] = {};
    Float_t l_onshell_phi[CONTROL] = {};
    Float_t l_onshell_mass[CONTROL] = {};

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

    Float_t fatjet_pt{0.0};
    Float_t fatjet_eta{0.0};
    Float_t fatjet_phi{0.0};
    Float_t fatjet_mass{0.0};

    Float_t Hbb_pt{0.0};
    Float_t Hbb_eta{0.0};
    Float_t Hbb_phi{0.0};
    Float_t Hbb_mass{0.0};

    Float_t bjet_resc_fact_1{0.0};
    Float_t bjet_resc_fact_2{0.0};
    Float_t fatjet_resc_fact{0.0};
    Float_t mh{0.0};
    Float_t mw{0.0};
    Float_t unclust_dpx{0.0};
    Float_t unclust_dpy{0.0};
    Float_t bjet_resc_dpx{0.0};
    Float_t bjet_resc_dpy{0.0};
    Float_t weight{0.0};
    Int_t num_sol{0};
};

#endif