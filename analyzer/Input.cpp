#include "Input.hpp"

EstimatorInput::EstimatorInput(std::vector<TLorentzVector>&& p4, std::string const& pdf_file)
: p4(std::move(p4))
{
    auto file_pdf = std::make_unique<TFile>(pdf_file.c_str(), "READ");

    size_t n_pdf1d = static_cast<size_t>(PDF1dSLRes::count);
    pdf1d.reserve(n_pdf1d);
    for (size_t i = 0; i < n_pdf1d; ++i)
    {
        auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get(pdf1d_names.at(i))));
        pdf1d.push_back(std::move(pdf));
    }

    size_t n_pdf2d = static_cast<size_t>(PDF2dSLRes::count);
    pdf2d.reserve(n_pdf2d);
    for (size_t i = 0; i < n_pdf2d; ++i)
    {
        auto pdf = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get(pdf2d_names.at(i))));
        pdf2d.push_back(std::move(pdf));
    }

    file_pdf->Close();
}