#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>
#include <map>

#include "TFile.h"
#include "TString.h"

#include "Event.hpp"
#include "Constants.hpp"
#include "Input.hpp"

class Analyzer
{
    private:
    std::map<TString, Channel> m_file_map;
    std::map<std::string, size_t> m_index;
    TString m_tree_name;
    TString m_pdf_file_name;

    public:
    Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name = "");
    void Analyze();
    void AnalyzeEvent(Event const& event, TTree* tree);

    private:
    HistVec1d_t PDFs1dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const;
    HistVec2d_t PDFs2dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const;
};


#endif