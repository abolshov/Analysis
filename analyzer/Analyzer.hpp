#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>

#include "TChain.h"
#include "TString.h"

#include "Event.hpp"

class Analyzer
{
    private:
    std::vector<TString> m_input_files;
    std::unique_ptr<TChain> m_chain;
    std::unique_ptr<Event> m_event;

    public:
    Analyzer(TString const& tree_name, std::vector<TString> const& input_files);
    void Analyze();
};


#endif