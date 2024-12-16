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

class Analyzer
{
    private:
    Storage m_storage;
    std::map<TString, Channel> m_file_map;
    TString m_tree_name;
    EstimatorSingLep m_esl; 

    public:
    Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode);
    void ProcessFile(TString const& name, Channel ch);
};


#endif