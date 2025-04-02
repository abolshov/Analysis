#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <vector>
#include <memory>
#include <map>
#include <fstream>

#include "TFile.h"
#include "TString.h"

#include "Constants.hpp"
#include "Definitions.hpp"
#include "Storage.hpp"
#include "Estimator.hpp"
#include "HistManager.hpp"
#include "Recorder.hpp"
#include "EstimatorData.hpp"

class Analyzer
{
    private:
    Storage m_storage;
    Estimator m_estimator;
    // HistManager m_hm; // eliminate ?
    Recorder m_recorder;
    std::unique_ptr<EstimatorData> m_estimator_data;
    bool m_record_iterations = false;
    bool m_record_output = true;

    inline static int counter = 0;

    UTree_t MakeTree();

    public:
    Analyzer(std::map<Channel, TString> const& pdf_file_map);
    
    void ProcessDataset(Dataset_t const& dataset); // new interface
    void ProcessSample(Sample const& sample); // new interface
    void ProcessEvent(ULong64_t evt, UTree_t& input_tree, UTree_t& output_tree, Channel ch, bool is_bkg); // new interface
    void ProcessEventSlice(ULong64_t evt, UTree_t& input_tree, UTree_t& output_tree); // estimate mass of the true combination if it exists
};


#endif