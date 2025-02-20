#include "EstimationRecorder.hpp"

EstimationRecorder::EstimationRecorder(TString const& out_file_name)
:   m_file(std::make_unique<TFile>(out_file_name, "RECREATE"))
{}