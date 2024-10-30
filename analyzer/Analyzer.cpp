#include "Analyzer.hpp"
#include "Definitions.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)     
{
    if (mode == Mode::Validation)
    {
        validator.Set1dPDF(PDFs1dFromFile(pdf_file_name, pdf1d_names));
        validator.Set2dPDF(PDFs2dFromFile(pdf_file_name, pdf2d_names));
    }
}

HistVec1d_t Analyzer::PDFs1dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec1d_t res;
    auto file = std::unique_ptr<TFile>(TFile::Open(file_name, "READ"));
    for (auto const& pdf_name: pdf_list)
    {
        std::unique_ptr<TH1F> hist(file->Get<TH1F>(pdf_name));
        hist->SetDirectory(nullptr); // unregister hist from root's garbage collector to avoid double free
        res.push_back(std::move(hist));
    }
    file->Close();
    return res;
}

HistVec2d_t Analyzer::PDFs2dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec2d_t res;
    auto file = std::unique_ptr<TFile>(TFile::Open(file_name, "READ"));
    for (auto const& pdf_name: pdf_list)
    {
        std::unique_ptr<TH2F> hist(file->Get<TH2F>(pdf_name));
        hist->SetDirectory(nullptr);
        res.push_back(std::move(hist));
    }
    file->Close();
    return res;
}


void Analyzer::Analyze()
{
    for (auto const& [file_name, channel]: m_file_map)
    {
        TFile* file = TFile::Open(file_name);
        TTree* tree = static_cast<TTree*>(file->Get(m_tree_name));

        TString channel_name = channel == Channel::SL ? "Single lepton" : "Double Lepton";
        m_obj_index = channel == Channel::SL ? GenTruthIdxMapSL : GenTruthIdxMapDL;
        std::cout << "File: " << file_name << ", channel: " << channel_name << "\n";

        m_pdf1d_index = PDF1dResolvedIdxMap;
        m_pdf2d_index = PDF2dResolvedIdxMap;

        Event event(tree, channel);
        // AnalyzeEvent(event, tree);

        validator.SetTruthIndex(m_obj_index);
        validator.SetPDF1dIndex(m_pdf1d_index);
        validator.SetPDF1dIndex(m_pdf2d_index);
        RunValidation(event, tree);

        file->Close();
    }
}

void Analyzer::AnalyzeEvent(Event const& event, TTree* tree)
{
    for (long long i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        size_t met_idx = m_obj_index.at("met");
        std::cout << "gen met: " << event.gen_truth.pt[met_idx] << " " << event.gen_truth.mass[met_idx] << "\n";
        std::cout << "reco jet: " << event.reco_jet.nRecoJet << "\n";
        std::cout << "gen jet: " << event.gen_jet.nGenJet << "\n";
        std::cout << "nu: " << event.nu.pdgId[0] << " " << event.nu.pdgId[1] << "\n";
        std::cout << "reco lep: " << event.reco_lep.lep_type[0] << " " << event.reco_lep.lep_iso[0] << " " << event.reco_lep.lep_iso[1] << "\n";
        std::cout << "lep1: " << event.reco_lep.pt[0] << " " << event.reco_lep.eta[0] << " " << event.reco_lep.phi[0] << " " << event.reco_lep.mass[0] << "\n";
        std::cout << "lep2: " << event.reco_lep.pt[1] << " " << event.reco_lep.eta[1] << " " << event.reco_lep.phi[1] << " " << event.reco_lep.mass[1] << "\n";
        std::cout << "\n";
    }
}


void Analyzer::RunValidation(Event const& event, TTree* tree)
{
    // run validation with 1d pdf
    bool use_2d_pdf = false;
    bool skip_failures = true;
    for (long long i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);
        validator.FillVariables(event);
        validator.Compare(i, use_2d_pdf, skip_failures);
    }

    validator.DrawStack({"raw_ratio_px_1", "corr_ratio_px_1"}, "Px ratio comarison for pair 1 with 1d PDF", "px_ratio_1_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_py_1", "corr_ratio_py_1"}, "Py ratio comarison for pair 1 with 1d PDF", "py_ratio_1_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_pz_1", "corr_ratio_pz_1"}, "Pz ratio comarison for pair 1 with 1d PDF", "pz_ratio_1_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_E_1", "corr_ratio_E_1"}, "E ratio comarison for pair 1 with 1d PDF", "E_ratio_1_pdf1d_cmp");

    validator.DrawStack({"raw_ratio_px_2", "corr_ratio_px_2"}, "Px ratio comarison for pair 2 with 1d PDF", "px_ratio_2_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_py_2", "corr_ratio_py_2"}, "Py ratio comarison for pair 2 with 1d PDF", "py_ratio_2_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_pz_2", "corr_ratio_pz_2"}, "Pz ratio comarison for pair 2 with 1d PDF", "pz_ratio_2_pdf1d_cmp");
    validator.DrawStack({"raw_ratio_E_2", "corr_ratio_E_2"}, "E ratio comarison for pair 2 with 1d PDF", "E_ratio_2_pdf1d_cmp");

    validator.DrawStack({"true_px_1", "corr_px_1", "raw_px_1"}, "Px comparison for pair 1 with 1d PDF", "px_1_pdf1d_cmp");
    validator.DrawStack({"true_py_1", "corr_py_1", "raw_py_1"}, "Py comparison for pair 1 with 1d PDF", "py_1_pdf1d_cmp");
    validator.DrawStack({"true_pz_1", "corr_pz_1", "raw_pz_1"}, "Pz comparison for pair 1 with 1d PDF", "pz_1_pdf1d_cmp");
    validator.DrawStack({"true_E_1", "corr_E_1", "raw_E_1"}, "E comparison for pair 1 with 1d PDF", "E_1_pdf1d_cmp");

    validator.DrawStack({"true_px_2", "corr_px_2", "raw_px_2"}, "Px comparison for pair 2 with 1d PDF", "px_2_pdf1d_cmp");
    validator.DrawStack({"true_py_2", "corr_py_2", "raw_py_2"}, "Py comparison for pair 2 with 1d PDF", "py_2_pdf1d_cmp");
    validator.DrawStack({"true_pz_2", "corr_pz_2", "raw_pz_2"}, "Pz comparison for pair 2 with 1d PDF", "pz_2_pdf1d_cmp");
    validator.DrawStack({"true_E_2", "corr_E_2", "raw_E_2"}, "E comparison for pair 2 with 1d PDF", "E_2_pdf1d_cmp");

    validator.DrawStack({"true_px_met", "corr_px_met", "raw_px_met"}, "Px comparison for MET with 1d PDF", "px_met_pdf1d_cmp");
    validator.DrawStack({"true_py_met", "corr_py_met", "raw_py_met"}, "Py comparison for MET with 1d PDF", "py_met_pdf1d_cmp");

    // reset histograms and run validation with 2d pdf
    validator.ResetHistograms();

    use_2d_pdf = true;
    for (long long i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);
        validator.FillVariables(event);
        validator.Compare(i, use_2d_pdf, skip_failures);
    }

    validator.DrawStack({"raw_ratio_px_1", "corr_ratio_px_1"}, "Px ratio comarison for pair 1 with 2d PDF", "px_ratio_1_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_py_1", "corr_ratio_py_1"}, "Py ratio comarison for pair 1 with 2d PDF", "py_ratio_1_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_pz_1", "corr_ratio_pz_1"}, "Pz ratio comarison for pair 1 with 2d PDF", "pz_ratio_1_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_E_1", "corr_ratio_E_1"}, "E ratio comarison for pair 1 with 2d PDF", "E_ratio_1_pdf2d_cmp");

    validator.DrawStack({"raw_ratio_px_2", "corr_ratio_px_2"}, "Px ratio comarison for pair 2 with 2d PDF", "px_ratio_2_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_py_2", "corr_ratio_py_2"}, "Py ratio comarison for pair 2 with 2d PDF", "py_ratio_2_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_pz_2", "corr_ratio_pz_2"}, "Pz ratio comarison for pair 2 with 2d PDF", "pz_ratio_2_pdf2d_cmp");
    validator.DrawStack({"raw_ratio_E_2", "corr_ratio_E_2"}, "E ratio comarison for pair 2 with 2d PDF", "E_ratio_2_pdf2d_cmp");

    validator.DrawStack({"true_px_1", "corr_px_1", "raw_px_1"}, "Px comparison for pair 1 with 2d PDF", "px_1_pdf2d_cmp");
    validator.DrawStack({"true_py_1", "corr_py_1", "raw_py_1"}, "Py comparison for pair 1 with 2d PDF", "py_1_pdf2d_cmp");
    validator.DrawStack({"true_pz_1", "corr_pz_1", "raw_pz_1"}, "Pz comparison for pair 1 with 2d PDF", "pz_1_pdf2d_cmp");
    validator.DrawStack({"true_E_1", "corr_E_1", "raw_E_1"}, "E comparison for pair 1 with 2d PDF", "E_1_pdf2d_cmp");

    validator.DrawStack({"true_px_2", "corr_px_2", "raw_px_2"}, "Px comparison for pair 2 with 2d PDF", "px_2_pdf2d_cmp");
    validator.DrawStack({"true_py_2", "corr_py_2", "raw_py_2"}, "Py comparison for pair 2 with 2d PDF", "py_2_pdf2d_cmp");
    validator.DrawStack({"true_pz_2", "corr_pz_2", "raw_pz_2"}, "Pz comparison for pair 2 with 2d PDF", "pz_2_pdf2d_cmp");
    validator.DrawStack({"true_E_2", "corr_E_2", "raw_E_2"}, "E comparison for pair 2 with 2d PDF", "E_2_pdf2d_cmp");

    validator.DrawStack({"true_px_met", "corr_px_met", "raw_px_met"}, "Px comparison for MET with 2d PDF", "px_met_pdf2d_cmp");
    validator.DrawStack({"true_py_met", "corr_py_met", "raw_py_met"}, "Py comparison for MET with 2d PDF", "py_met_pdf2d_cmp");
}