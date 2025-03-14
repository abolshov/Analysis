#include "EstimationRecorder.hpp"

EstimationRecorder::EstimationRecorder(TString const& dbg_file_name)
{
    if (!dbg_file_name.EqualTo(""))
    {
        m_file = std::make_unique<TFile>(dbg_file_name, "RECREATE");
    }
}

void EstimationRecorder::OpenFile(TString const& name)
{
    if (m_file != nullptr)
    {
        m_file->Close();
    }
    m_file = std::make_unique<TFile>(name, "RECREATE");
}
