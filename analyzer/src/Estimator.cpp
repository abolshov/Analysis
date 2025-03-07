#include "Estimator.hpp"

Estimator::Estimator(TString const& pdf_file_name_sl, TString const& pdf_file_name_dl, TString const& dbg_file_name)
:   m_estimator_sl(pdf_file_name_sl, dbg_file_name)
,   m_estimator_dl(pdf_file_name_dl, dbg_file_name)
{}

std::optional<Float_t> Estimator::EstimateMass(VecLVF_t const& jets, VecLVF_t const& leptons, LorentzVectorF_t const& met, ULong64_t evt_id, Channel ch)
{
    if (ch == Channel::SL)
    {
        return m_estimator_sl.EstimateMass(jets, leptons, met, evt_id);
    }
    else if (ch == Channel::DL)
    {
        return m_estimator_dl.EstimateMass(jets, leptons, met, evt_id);
    }
    else 
    {
        throw std::runtime_error("Attempting to process data in unknown channel");
    }
}
