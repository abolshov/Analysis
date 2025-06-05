#include "EvaluatorBase.hpp"

#include <vector>

EvaluatorBase::EvaluatorBase(std::ifstream& file_stream, Channel ch)
:   m_chain(std::make_unique<TChain>("Events"))
,   m_buf(std::make_unique<Buffer>())
,   m_out_file(std::make_unique<TFile>(ch == Channel::DL ? "pdfs_dl.root" : "pdfs_sl.root", "RECREATE"))
{
    std::string file_name;
    std::vector<std::string> file_names;
    while (file_stream >> file_name)
    {
        file_names.emplace_back(std::move(file_name));
        for (auto const& file_name: file_names)
        {
            m_chain->Add(file_name.c_str());
        }
    }
    m_buf->ConnectTree(m_chain.get(), ch);
}
