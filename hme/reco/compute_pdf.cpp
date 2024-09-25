#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TString.h"

inline constexpr size_t N_RECO_JETS = 20;
inline constexpr size_t N_BINS = 2000;

template <class T>
double GetPDFScaleFactor(std::unique_ptr<T> const& hist)
{
    int binmax = hist->GetMaximumBin();
    return hist->GetBinContent(binmax);
}

int main()
{
    std::vector<TString> input_files;
    std::ifstream files("files.txt");
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

    chain->SetBranchAddress("ncentralJet", &ncentralJet);
    chain->SetBranchAddress("centralJet_pt", centralJet_pt);
    chain->SetBranchAddress("centralJet_eta", centralJet_eta);
    chain->SetBranchAddress("centralJet_phi", centralJet_phi);
    chain->SetBranchAddress("centralJet_mass", centralJet_mass);

    // V2: hadronic
    Double_t        genV2prod1_pt;
    Double_t        genV2prod1_eta;
    Double_t        genV2prod1_phi;
    Double_t        genV2prod1_mass;

    Double_t        genV2prod2_pt;
    Double_t        genV2prod2_eta;
    Double_t        genV2prod2_phi;
    Double_t        genV2prod2_mass;

    // H->bb
    Double_t        genb1_pt;
    Double_t        genb1_eta;
    Double_t        genb1_phi;
    Double_t        genb1_mass;

    Double_t        genb2_pt;
    Double_t        genb2_eta;
    Double_t        genb2_phi;
    Double_t        genb2_mass;

    chain->SetBranchAddress("genV2prod1_pt", &genV2prod1_pt);
    chain->SetBranchAddress("genV2prod1_eta", &genV2prod1_eta);
    chain->SetBranchAddress("genV2prod1_phi", &genV2prod1_phi);
    chain->SetBranchAddress("genV2prod1_mass", &genV2prod1_mass);

    chain->SetBranchAddress("genV2prod2_pt", &genV2prod2_pt);
    chain->SetBranchAddress("genV2prod2_eta", &genV2prod2_eta);
    chain->SetBranchAddress("genV2prod2_phi", &genV2prod2_phi);
    chain->SetBranchAddress("genV2prod2_mass", &genV2prod2_mass);

    chain->SetBranchAddress("genb1_pt", &genb1_pt);
    chain->SetBranchAddress("genb1_eta", &genb1_eta);
    chain->SetBranchAddress("genb1_phi", &genb1_phi);
    chain->SetBranchAddress("genb1_mass", &genb1_mass);
    chain->SetBranchAddress("genb2_pt", &genb2_pt);
    chain->SetBranchAddress("genb2_eta", &genb2_eta);
    chain->SetBranchAddress("genb2_phi", &genb2_phi);
    chain->SetBranchAddress("genb2_mass", &genb2_mass);

    auto pdf_b1 = std::make_unique<TH1F>("pdf_b1", "1d PDF for leading b jet correction", N_BINS, 0, 8);
    auto pdf_b2 = std::make_unique<TH1F>("pdf_b2", "1d PDF for subleading b jet correction", N_BINS, 0, 8);
    auto pdf_b1b2 = std::make_unique<TH2F>("pdf_b1b2", "2d PDF simultaneous b jet corrrections", N_BINS, 0, 8, N_BINS, 0, 8);

    auto pdf_mbb = std::make_unique<TH1F>("pdf_mbb", "1d PDF of H->bb mass with true corrections applied", 5*N_BINS, 0, 200);

    size_t nEvents = chain->GetEntries();
    for (size_t i = 0; i < nEvents; ++i)
    {
        chain->GetEntry(i);

        if (ncentralJet < 4)
        {
            continue;
        }

        TLorentzVector genb1_p4, genb2_p4, genq1_p4, genq2_p4; // quarks
        genb1_p4.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        genb2_p4.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);
        genq1_p4.SetPtEtaPhiM(genV2prod1_pt, genV2prod1_eta, genV2prod1_phi, genV2prod1_mass);
        genq2_p4.SetPtEtaPhiM(genV2prod2_pt, genV2prod2_eta, genV2prod2_phi, genV2prod2_mass);

        if (genb1_p4.DeltaR(genb2_p4) < 0.4 || genq1_p4.DeltaR(genq2_p4) < 0.4)
        {
            continue;
        }
        
        // establish jet with greater pt for computing corrections
        // jets are sorted by btag so first two are b jets
        TLorentzVector reco_bj1_p4, reco_bj2_p4;
        if (centralJet_pt[0] > centralJet_pt[1])
        {
            reco_bj1_p4.SetPtEtaPhiM(centralJet_pt[0], centralJet_eta[0], centralJet_phi[0], centralJet_mass[0]);
            reco_bj2_p4.SetPtEtaPhiM(centralJet_pt[1], centralJet_eta[1], centralJet_phi[1], centralJet_mass[1]);
        }
        else 
        {
            reco_bj2_p4.SetPtEtaPhiM(centralJet_pt[0], centralJet_eta[0], centralJet_phi[0], centralJet_mass[0]);
            reco_bj1_p4.SetPtEtaPhiM(centralJet_pt[1], centralJet_eta[1], centralJet_phi[1], centralJet_mass[1]);
        }

        TLorentzVector lead_quark = genb1_p4.Pt() > genb2_p4.Pt() ? genb1_p4 : genb2_p4;
        TLorentzVector subl_quark = genb1_p4.Pt() > genb2_p4.Pt() ? genb2_p4 : genb1_p4;
        double c1 = lead_quark.Pt()/reco_bj1_p4.Pt();
        double c2 = subl_quark.Pt()/reco_bj2_p4.Pt();
        pdf_b1->Fill(c1);
        pdf_b2->Fill(c1);
        pdf_b1b2->Fill(c1, c2); 
        pdf_mbb->Fill((c1*reco_bj1_p4 + c2*reco_bj2_p4).M());
    }

    pdf_b1->Scale(1.0/GetPDFScaleFactor(pdf_b1));
    pdf_b2->Scale(1.0/GetPDFScaleFactor(pdf_b2));
    pdf_mbb->Scale(1.0/GetPDFScaleFactor(pdf_mbb));

    auto output = std::make_unique<TFile>("pdf.root","RECREATE");
    pdf_b1->Write();
    pdf_b2->Write();
    pdf_mbb->Write();
    pdf_b1b2->Write();
	output->Write();
	output->Close();

    return 0;
}