#include "Estimator.hpp"

Estimator::Estimator(TString const& pdf_file_name_sl, TString const& pdf_file_name_dl, AggregationMode aggr_mode)
:   m_estimator_sl(pdf_file_name_sl, aggr_mode)
,   m_estimator_dl(pdf_file_name_dl, aggr_mode)
{}

OptArrF_t<ESTIM_OUT_SZ> Estimator::EstimateMass(VecLVF_t const& jets, std::vector<Float_t> const& resolutions, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id, Channel ch)
{
    if (ch == Channel::SL)
    {
        return m_estimator_sl.EstimateMass(jets, resolutions, leptons, met, evt_id);
    }
    else if (ch == Channel::DL)
    {
        return m_estimator_dl.EstimateMass(jets, resolutions, leptons, met, evt_id);
    }
    else 
    {
        throw std::runtime_error("Attempting to process data in unknown channel");
    }
}

OptArrF_t<ESTIM_OUT_SZ> Estimator::EstimateMass(Event const& event, Channel ch)
{
    if (ch == Channel::SL)
    {
        return m_estimator_sl.EstimateMass(event);
    }
    else if (ch == Channel::DL)
    {
        return m_estimator_dl.EstimateMass(event);
    }
    else 
    {
        throw std::runtime_error("Attempting to process data in unknown channel");
    }
}

void Estimator::OpenDbgFile(TString const& dbg_file_name, Channel ch)
{
    if (ch == Channel::SL)
    {
        m_estimator_sl.OpenDbgFile(dbg_file_name);
    }
    else if (ch == Channel::DL)
    {
        m_estimator_dl.OpenDbgFile(dbg_file_name);
    }
    else 
    {
        throw std::runtime_error("Attempting to process data in unknown channel");
    }
}
