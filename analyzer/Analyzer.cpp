#include "Analyzer.hpp"
#include "Definitions.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)  
,   m_esl(pdf_file_name)   
{
    if (mode == Mode::Validation)
    {
        std::cout << "Constructing analyzer in validation mode\n";
    }
}