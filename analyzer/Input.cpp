#include "Input.hpp"

EstimatorInput::EstimatorInput(std::vector<TLorentzVector>&& p4, std::string const& pdf_file)
: p4(std::move(p4)),
  pdf1d(PDF1dSLRes::count),
  pdf2d(PDF2dSLRes::count)
{
    auto file_pdf = std::make_unique<TFile>(pdf_file.c_str(), "READ");

    for (int i = 0; i < PDF1dSLRes::count; ++i)
    {
        auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get(pdf1d_names.at(i))));
        pdf1d.push_back(std::move(pdf));
    }

    for (int i = 0; i < PDF2dSLRes::count; ++i)
    {
        auto pdf = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get(pdf2d_names.at(i))));
        pdf2d.push_back(std::move(pdf));
    }

    file_pdf.Close();
}