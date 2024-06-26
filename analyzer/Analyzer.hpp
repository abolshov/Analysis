#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>

#include "TChain.h"
#include "TString.h"

#include "EventData.hpp"
#include "Selector.hpp"

class Analyzer
{
    private:
        std::unique_ptr<TChain> m_chain;
        std::unique_ptr<EventData> m_data;
        std::vector<TString> m_input_files;
        Selector m_selector;

    public:
        Analyzer(TString const& tree_name, std::vector<TString> const& input_files);
        void Analyze();
};


#endif