#include "Recorder.hpp"

Recorder::Recorder(TString const& file_name)
{
    if (!file_name.EqualTo(""))
    {
        m_file = std::make_unique<TFile>(file_name, "RECREATE");
    }
}

void Recorder::OpenFile(TString const& file_name)
{
    if (m_file != nullptr)
    {
        m_file->Close();
    }
    m_file = std::make_unique<TFile>(file_name, "RECREATE");
}
