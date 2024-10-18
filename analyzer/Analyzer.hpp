#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>
#include <map>

#include "TChain.h"
#include "TString.h"

#include "Event.hpp"
#include "Constants.hpp"

class Analyzer
{
    private:
    std::map<TString, Channel> m_file_map;
    std::map<std::string, size_t> m_index;
    TString m_tree_name;

    public:
    Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map);
    void Analyze();
    void AnalyzeEvent(Event const& event, TTree* tree);
};


#endif