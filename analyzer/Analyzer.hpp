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

class Analyzer
{
    private:
    Storage m_storage;
    std::map<TString, Channel> m_file_map;
    TString m_tree_name;

    public:
    Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode);
    void ProcessFile(TString const& name, Channel ch);

    private:
    HistVec1d_t PDFs1dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const;
    HistVec2d_t PDFs2dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const;

    VecLVF_t GetRecoJetP4();
    VecLVF_t GetRecoLepP4(Channel ch);
    std::vector<Float_t> GetPNetRes();

    inline LorentzVectorF_t GetRecoMET() { return LorentzVectorF_t(m_storage.reco_met_pt, 0.0, m_storage.reco_met_phi, 0.0); }
    inline LorentzVectorF_t GetGenMET() { return LorentzVectorF_t(m_storage.gen_met_pt, 0.0, m_storage.gen_met_phi, 0.0); }

    VecLVF_t GetGenJetP4();
    VecLVF_t GetGenLepP4();
};


#endif