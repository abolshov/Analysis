#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>
#include <map>

#include "TFile.h"
#include "TString.h"

#include "Constants.hpp"
#include "Definitions.hpp"
#include "Storage.hpp"
#include "Estimator.hpp"
#include "HistManager.hpp"

class Analyzer
{
    private:
    Storage m_storage;
    std::map<TString, Channel> m_file_map;
    TString m_tree_name;
    EstimatorSingLep m_estimator; 
    HistManager m_hm;

    public:
    Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode);
    
    void ProcessFile(TString const& name, Channel ch);
    void ProcessEvent(ULong64_t evt, TTree* tree, Channel ch);
};


#endif