#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include "EstimatorSingleLep.hpp"
#include "EstimatorDoubleLep.hpp"

class Estimator
{
    public:
    Estimator(TString const& pdf_file_name_sl, TString const& pdf_file_name_dl, AggregationMode aggr_mode);
    OptArrF_t<ESTIM_OUT_SZ> EstimateMass(VecLVF_t const& jets, std::vector<Float_t> const& resolutions, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id, Channel ch);
    OptArrF_t<ESTIM_OUT_SZ> EstimateMass(Event const& event, Channel ch);
    void OpenDbgFile(TString const& dbg_file_name, Channel ch);
    inline EstimatorBase& GetEstimator(Channel ch) { return ch == Channel::SL ? static_cast<EstimatorBase&>(m_estimator_sl) : static_cast<EstimatorBase&>(m_estimator_dl); }

    private:
    EstimatorSingleLep m_estimator_sl;
    EstimatorDoubleLep m_estimator_dl;
};

#endif