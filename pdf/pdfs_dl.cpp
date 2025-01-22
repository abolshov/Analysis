#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <fstream>
#include <unordered_set>
#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"

inline constexpr size_t N_RECO_JETS = 20;
inline constexpr size_t N_BINS = 100;

template <class T>
double GetPDFScaleFactor(std::unique_ptr<T> const& hist)
{
    int binmax = hist->GetMaximumBin();
    return hist->GetBinContent(binmax);
}

int FindBestMatch(TLorentzVector const& quark, std::vector<TLorentzVector> const& jets)
{
    int res = -1;
    int sz = jets.size();
    double min_dr = 10.0;
    for (int i = 0; i < sz; ++i)
    {
        double dr = jets[i].DeltaR(quark);
        if (dr < min_dr)
        {
            min_dr = dr;
            res = i;
        }
    }
    if (min_dr < 0.4)
    {
        return res;
    }
    return -1;
}

bool CorrectLepReco(int lep_type, int lep_genLep_kind)
{
    bool reco_lep_mu = (lep_type == 2);
    bool reco_lep_ele = (lep_type == 1);

    bool gen_lep_mu = ((lep_genLep_kind == 2) || (lep_genLep_kind == 4));
    bool gen_lep_ele = ((lep_genLep_kind == 1) || (lep_genLep_kind == 3));

    bool corr_lep_reco = ((reco_lep_mu && gen_lep_mu) || (reco_lep_ele && gen_lep_ele));
    return corr_lep_reco;
}

int main()
{
    std::vector<TString> input_files;
    std::ifstream files("files_dl.txt");
    std::string fname;
    while (files >> fname)
    {
        input_files.push_back(std::move(fname));
    }

    TString tree_name = "Events";
    auto chain = std::make_unique<TChain>(tree_name);
    
    for (auto const& input_file: input_files)
    {
        std::cout << "Adding " << input_file << "\n";
        chain->Add(input_file);
    }

    Int_t           ncentralJet;
    Float_t         centralJet_pt[N_RECO_JETS];   
    Float_t         centralJet_eta[N_RECO_JETS];   
    Float_t         centralJet_phi[N_RECO_JETS];   
    Float_t         centralJet_mass[N_RECO_JETS]; 
    Float_t         centralJet_PNetRegPtRawCorr[N_RECO_JETS];    
    Float_t         centralJet_PNetRegPtRawRes[N_RECO_JETS]; 

    chain->SetBranchAddress("ncentralJet", &ncentralJet);
    chain->SetBranchAddress("centralJet_pt", centralJet_pt);
    chain->SetBranchAddress("centralJet_eta", centralJet_eta);
    chain->SetBranchAddress("centralJet_phi", centralJet_phi);
    chain->SetBranchAddress("centralJet_mass", centralJet_mass);
    chain->SetBranchAddress("centralJet_PNetRegPtRawCorr", centralJet_PNetRegPtRawCorr);
    chain->SetBranchAddress("centralJet_PNetRegPtRawRes", centralJet_PNetRegPtRawRes);

    Float_t         lep1_pt;
    Float_t         lep1_eta;
    Float_t         lep1_phi;
    Float_t         lep1_mass;

    Int_t           lep1_type;
    Int_t           lep1_genLep_kind;

    Float_t         lep2_pt;
    Float_t         lep2_eta;
    Float_t         lep2_phi;
    Float_t         lep2_mass;

    Int_t           lep2_type;
    Int_t           lep2_genLep_kind;

    chain->SetBranchAddress("lep1_pt", &lep1_pt);
    chain->SetBranchAddress("lep1_eta", &lep1_eta);
    chain->SetBranchAddress("lep1_phi", &lep1_phi);
    chain->SetBranchAddress("lep1_mass", &lep1_mass);

    chain->SetBranchAddress("lep1_type", &lep1_type);
    chain->SetBranchAddress("lep1_genLep_kind", &lep1_genLep_kind);

    chain->SetBranchAddress("lep2_pt", &lep2_pt);
    chain->SetBranchAddress("lep2_eta", &lep2_eta);
    chain->SetBranchAddress("lep2_phi", &lep2_phi);
    chain->SetBranchAddress("lep2_mass", &lep2_mass);

    chain->SetBranchAddress("lep2_type", &lep2_type);
    chain->SetBranchAddress("lep2_genLep_kind", &lep2_genLep_kind);

    Float_t        genHVV_pt;
    Float_t        genHVV_eta;
    Float_t        genHVV_phi;
    Float_t        genHVV_mass;    

    Float_t        genHbb_pt;
    Float_t        genHbb_eta;
    Float_t        genHbb_phi;
    Float_t        genHbb_mass;

    chain->SetBranchAddress("genHVV_pt", &genHVV_pt);
    chain->SetBranchAddress("genHVV_eta", &genHVV_eta);
    chain->SetBranchAddress("genHVV_phi", &genHVV_phi);
    chain->SetBranchAddress("genHVV_mass", &genHVV_mass);

    chain->SetBranchAddress("genHbb_pt", &genHbb_pt);
    chain->SetBranchAddress("genHbb_eta", &genHbb_eta);
    chain->SetBranchAddress("genHbb_phi", &genHbb_phi);
    chain->SetBranchAddress("genHbb_mass", &genHbb_mass);    

    // V2: hadronic
    Float_t        genV2_mass;
    Float_t        genV1_mass;

    Float_t        genV2prod1_pt;
    Float_t        genV2prod1_eta;
    Float_t        genV2prod1_phi;
    Float_t        genV2prod1_mass;

    Float_t        genV2prod2_pt;
    Float_t        genV2prod2_eta;
    Float_t        genV2prod2_phi;
    Float_t        genV2prod2_mass;

    // H->bb
    Float_t        genb1_pt;
    Float_t        genb1_eta;
    Float_t        genb1_phi;
    Float_t        genb1_mass;

    Float_t        genb2_pt;
    Float_t        genb2_eta;
    Float_t        genb2_phi;
    Float_t        genb2_mass;

    Float_t         PuppiMET_pt;
    Float_t         PuppiMET_phi;

    // neutrino
    Float_t        genV1prod2_pt;
    Float_t        genV1prod2_eta;
    Float_t        genV1prod2_phi;
    Float_t        genV1prod2_mass;

    chain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    chain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi);

    chain->SetBranchAddress("genV2prod1_pt", &genV2prod1_pt);
    chain->SetBranchAddress("genV2prod1_eta", &genV2prod1_eta);
    chain->SetBranchAddress("genV2prod1_phi", &genV2prod1_phi);
    chain->SetBranchAddress("genV2prod1_mass", &genV2prod1_mass);

    chain->SetBranchAddress("genV2prod2_pt", &genV2prod2_pt);
    chain->SetBranchAddress("genV2prod2_eta", &genV2prod2_eta);
    chain->SetBranchAddress("genV2prod2_phi", &genV2prod2_phi);
    chain->SetBranchAddress("genV2prod2_mass", &genV2prod2_mass);

    chain->SetBranchAddress("genV1prod2_pt", &genV1prod2_pt);
    chain->SetBranchAddress("genV1prod2_eta", &genV1prod2_eta);
    chain->SetBranchAddress("genV1prod2_phi", &genV1prod2_phi);
    chain->SetBranchAddress("genV1prod2_mass", &genV1prod2_mass);

    chain->SetBranchAddress("genb1_pt", &genb1_pt);
    chain->SetBranchAddress("genb1_eta", &genb1_eta);
    chain->SetBranchAddress("genb1_phi", &genb1_phi);
    chain->SetBranchAddress("genb1_mass", &genb1_mass);
    chain->SetBranchAddress("genb2_pt", &genb2_pt);
    chain->SetBranchAddress("genb2_eta", &genb2_eta);
    chain->SetBranchAddress("genb2_phi", &genb2_phi);
    chain->SetBranchAddress("genb2_mass", &genb2_mass);

    chain->SetBranchAddress("genV2_mass", &genV2_mass);
    chain->SetBranchAddress("genV1_mass", &genV1_mass);

    auto pdf_b1b2 = std::make_unique<TH2F>("pdf_b1b2", "2d PDF simultaneous b jet corrrections", N_BINS, 0, 8, N_BINS, 0, 8);
    auto pdf_mbb = std::make_unique<TH1F>("pdf_mbb", "1d PDF of H->bb mass with true corrections applied", N_BINS, 0, 200);
    auto pdf_b1 = std::make_unique<TH1F>("pdf_b1", "1d PDF for leading b jet correction", N_BINS, 0, 8);
    auto pdf_b2 = std::make_unique<TH1F>("pdf_b2", "1d PDF for subleading b jet correction", N_BINS, 0, 8);
    auto pdf_mw_onshell = std::make_unique<TH1F>("pdf_mw_onshell", "1d PDF of onshell W", N_BINS, 0, 100);
    auto pdf_mw_offshell = std::make_unique<TH1F>("pdf_mw_offshell", "1d PDF of offshell W", N_BINS, 0, 100);
    
    long long nEvents = chain->GetEntries();
    for (long long i = 0; i < nEvents; ++i)
    {
        chain->GetEntry(i);

        if (ncentralJet < 2)
        {
            continue;
        }

        TLorentzVector genb1_p4, genb2_p4;
        genb1_p4.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        genb2_p4.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);

        if (genb1_p4.DeltaR(genb2_p4) < 0.4)
        {
            continue;
        }

        std::vector<TLorentzVector> jets;
        for (int j = 0; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);
            jets.push_back(jet);
        }

        bool correct_leptons = CorrectLepReco(lep1_type, lep1_genLep_kind) && CorrectLepReco(lep2_type, lep2_genLep_kind);
        if (!correct_leptons)
        {
            continue;
        }

        std::unordered_set<int> match_idx;
        int b1_match = FindBestMatch(genb1_p4, jets);
        int b2_match = FindBestMatch(genb2_p4, jets);
        match_idx.insert(b1_match);
        match_idx.insert(b2_match);

        if (match_idx.size() != 2 || match_idx.count(-1))
        {
            continue;
        }

        TLorentzVector const& reco_bj1_p4 = jets[b1_match];
        TLorentzVector const& reco_bj2_p4 = jets[b2_match];

        double c1 = genb1_p4.Pt()/reco_bj1_p4.Pt();
        double c2 = genb2_p4.Pt()/reco_bj2_p4.Pt();

        if (genV1_mass > genV2_mass)
        {
            pdf_mw_onshell->Fill(genV1_mass);
            pdf_mw_offshell->Fill(genV2_mass);
        }
        else 
        {
            pdf_mw_onshell->Fill(genV2_mass);
            pdf_mw_offshell->Fill(genV1_mass);
        }

        pdf_b1->Fill(c1);
        pdf_b2->Fill(c2);
        pdf_b1b2->Fill(c1, c2); 
        pdf_mbb->Fill((c1*reco_bj1_p4 + c2*reco_bj2_p4).M());
    }

    std::vector<double> xval = {0.0,  0.02, 0.04, 0.06, 0.08, 0.1,  0.12, 0.14, 0.16, 0.18, 0.2,  0.22, 0.24, 0.26, 0.28, 0.3,  0.32, 0.34, 0.36,
        0.38, 0.4,  0.42, 0.44, 0.46, 0.48, 0.5,  0.52, 0.54, 0.56, 0.58, 0.6,  0.62, 0.64, 0.66, 0.68, 0.7,  0.72, 0.74,
        0.76, 0.78, 0.8,  0.82, 0.84, 0.86, 0.88, 0.9,  0.92, 0.94, 0.96, 0.98, 1.0,  1.02, 1.04, 1.06, 1.08, 1.1,  1.12,
        1.14, 1.16, 1.18, 1.2,  1.22, 1.24, 1.26, 1.28, 1.3,  1.32, 1.34, 1.36, 1.38, 1.4,  1.42, 1.44, 1.46, 1.48, 1.5,
        1.52, 1.54, 1.56, 1.58, 1.6,  1.62, 1.64, 1.66, 1.68, 1.7,  1.72, 1.74, 1.76, 1.78, 1.8,  1.82, 1.84, 1.86, 1.88,
        1.9,  1.92, 1.94, 1.96, 1.98, 2.0,  2.02, 2.04, 2.06, 2.08, 2.1,  2.12, 2.14, 2.16, 2.18, 2.2,  2.22, 2.24, 2.26,
        2.28, 2.3,  2.32, 2.34, 2.36, 2.38, 2.4,  2.42, 2.44, 2.46, 2.48, 2.5,  2.52, 2.54, 2.56, 2.58, 2.6,  2.62, 2.64,
        2.66, 2.68, 2.7,  2.72, 2.74, 2.76, 2.78, 2.8,  2.82, 2.84, 2.86, 2.88, 2.9,  2.92, 2.94, 2.96, 2.98, 3.0,  3.02,
        3.04, 3.06, 3.08, 3.1,  3.12, 3.14, 3.16, 3.18, 3.2,  3.22, 3.24, 3.26, 3.28, 3.3,  3.32, 3.34, 3.36, 3.38, 3.4,
        3.42, 3.44, 3.46, 3.48, 3.5,  3.52, 3.54, 3.56, 3.58, 3.6,  3.62, 3.64, 3.66, 3.68, 3.7,  3.72, 3.74, 3.76, 3.78,
        3.8,  3.82, 3.84, 3.86, 3.88, 3.9,  3.92, 3.94, 3.96, 3.98, 4.0,  4.02, 4.04, 4.06, 4.08, 4.1,  4.12, 4.14, 4.16,
        4.18, 4.2,  4.22, 4.24, 4.26, 4.28, 4.3,  4.32, 4.34, 4.36, 4.38, 4.4,  4.42, 4.44, 4.46, 4.48, 4.5,  4.52, 4.54,
        4.56, 4.58, 4.6,  4.62, 4.64, 4.66, 4.68, 4.7,  4.72, 4.74, 4.76, 4.78, 4.8,  4.82, 4.84, 4.86, 4.88, 4.9,  4.92,
        4.94, 4.96, 4.98, 5.0,  5.02, 5.04, 5.06, 5.08, 5.1,  5.12, 5.14, 5.16, 5.18, 5.2,  5.22, 5.24, 5.26, 5.28, 5.3,
        5.32, 5.34, 5.36, 5.38, 5.4,  5.42, 5.44, 5.46, 5.48, 5.5,  5.52, 5.54, 5.56, 5.58, 5.6,  5.62, 5.64, 5.66, 5.68,
        5.7,  5.72, 5.74, 5.76, 5.78, 5.8,  5.82, 5.84, 5.86, 5.88, 5.9,  5.92, 5.94, 5.96, 5.98, 6.0};

    std::vector<double> yval = {0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   2.0,   1.0,   4.0,   6.0,   7.0,   4.0,   4.0,   4.0,   9.0,
        6.0,   16.0,  8.0,   7.0,   8.0,   5.0,   6.0,   5.0,   4.0,   8.0,   14.0,  7.0,   21.0,  9.0,   7.0,   14.0,
        15.0,  16.0,  9.0,   19.0,  17.0,  28.0,  24.0,  40.0,  51.0,  58.0,  73.0,  88.0,  126.0, 173.0, 269.0, 371.0,
        474.0, 594.0, 695.0, 702.0, 777.0, 735.0, 742.0, 636.0, 593.0, 467.0, 458.0, 392.0, 383.0, 341.0, 319.0, 293.0,
        270.0, 239.0, 204.0, 184.0, 154.0, 151.0, 153.0, 133.0, 127.0, 101.0, 104.0, 120.0, 77.0,  70.0,  61.0,  57.0,
        74.0,  57.0,  73.0,  59.0,  56.0,  47.0,  30.0,  24.0,  38.0,  46.0,  33.0,  32.0,  21.0,  29.0,  30.0,  21.0,
        18.0,  25.0,  20.0,  17.0,  19.0,  6.0,   11.0,  14.0,  14.0,  9.0,   12.0,  4.0,   10.0,  11.0,  7.0,   5.0,
        7.0,   4.0,   5.0,   4.0,   8.0,   3.0,   2.0,   0.0,   2.0,   8.0,   6.0,   5.0,   0.0,   2.0,   2.0,   6.0,
        2.0,   1.0,   1.0,   1.0,   0.0,   2.0,   4.0,   0.0,   1.0,   2.0,   0.0,   2.0,   1.0,   2.0,   3.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   1.0,   0.0,   0.0,   0.0,   1.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0};

    auto pdf_b1_run2 = std::make_unique<TH1F>("pdf_b1_run2", "1d PDF for leading b jet correction", xval.size(), 0, xval.back());
    for (size_t i = 0; i < xval.size(); ++i)
    {
        pdf_b1_run2->SetBinContent(i+1, yval[i]);
    }

    pdf_b1->Scale(1.0/GetPDFScaleFactor(pdf_b1));
    pdf_b2->Scale(1.0/GetPDFScaleFactor(pdf_b2));
    pdf_mbb->Scale(1.0/GetPDFScaleFactor(pdf_mbb));
    pdf_b1b2->Scale(1.0/GetPDFScaleFactor(pdf_b1b2));
    pdf_mw_onshell->Scale(1.0/GetPDFScaleFactor(pdf_mw_onshell));
    pdf_mw_offshell->Scale(1.0/GetPDFScaleFactor(pdf_mw_offshell));

    auto output = std::make_unique<TFile>("pdf_dl.root", "RECREATE");
    pdf_b1->Write();
    pdf_b2->Write();
    pdf_mbb->Write();
    pdf_b1b2->Write();
    pdf_mw_onshell->Write();
    pdf_mw_offshell->Write();
    pdf_b1_run2->Write();
	output->Write();
	output->Close();

    return 0;
}