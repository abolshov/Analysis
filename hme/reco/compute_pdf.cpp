#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TChain.h"
#include "TString.h"

template <typename T>
std::vector<int> sort_indices(T* v, int n) 
{
  // initialize original index locations
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2];});

  return idx;
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
    Float_t         centralJet_pt[20];   //[ncentralJet]
    Float_t         centralJet_eta[20];   //[ncentralJet]
    Float_t         centralJet_phi[20];   //[ncentralJet]
    Float_t         centralJet_mass[20];   //[ncentralJet]

    Float_t         centralJet_btagPNetB[20];   //[ncentralJet]
    Float_t         centralJet_btagDeepFlavB[20];   //[ncentralJet]

    Double_t        genb1_pt;
    Double_t        genb1_eta;
    Double_t        genb1_phi;
    Double_t        genb1_mass;

    Double_t        genb2_pt;
    Double_t        genb2_eta;
    Double_t        genb2_phi;
    Double_t        genb2_mass;

    chain->SetBranchAddress("ncentralJet", &ncentralJet);
    chain->SetBranchAddress("centralJet_pt", &centralJet_pt);
    chain->SetBranchAddress("centralJet_eta", &centralJet_eta);
    chain->SetBranchAddress("centralJet_phi", &centralJet_phi);
    chain->SetBranchAddress("centralJet_mass", &centralJet_mass);

    chain->SetBranchAddress("centralJet_btagPNetB", &centralJet_btagPNetB);
    chain->SetBranchAddress("centralJet_btagDeepFlavB", &centralJet_btagDeepFlavB);

    chain->SetBranchAddress("genb1_pt", &genb1_pt);
    chain->SetBranchAddress("genb1_eta", &genb1_eta);
    chain->SetBranchAddress("genb1_phi", &genb1_phi);
    chain->SetBranchAddress("genb1_mass", &genb1_mass);

    chain->SetBranchAddress("genb2_pt", &genb2_pt);
    chain->SetBranchAddress("genb2_eta", &genb2_eta);
    chain->SetBranchAddress("genb2_phi", &genb2_phi);
    chain->SetBranchAddress("genb2_mass", &genb2_mass);

    auto pdf = std::make_unique<TH1F>("pdf", "pdf", 100, 0, 6);

    int nEvents = chain->GetEntries();
    for (int i = 0; i < nEvents; ++i)
    {
        chain->GetEntry(i);
        std::vector<int> sorted_by_PNet = sort_indices<Float_t>(centralJet_btagPNetB, ncentralJet);
        // std::vector<int> sorted_by_DeepFlav = sort_indices<Float_t>(centralJet_btagDeepFlavB, ncentralJet);

        int j1 = sorted_by_PNet[0];
        int j2 = sorted_by_PNet[1];

        TLorentzVector bj1, bj2;
        bj1.SetPtEtaPhiM(centralJet_pt[j1], centralJet_eta[j1], centralJet_phi[j1], centralJet_mass[j1]);
        bj2.SetPtEtaPhiM(centralJet_pt[j2], centralJet_eta[j2], centralJet_phi[j2], centralJet_mass[j2]);

        TLorentzVector bq1, bq2;
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);

        bool first_match = bq1.DeltaR(bj1) < 0.4 || bq1.DeltaR(bj2) < 0.4;
        bool second_match = bq2.DeltaR(bj1) < 0.4 || bq2.DeltaR(bj2) < 0.4;

        if (first_match && second_match)
        {
            TLorentzVector lead_jet = bj1.Pt() > bj2.Pt() ? bj1 : bj2;
            TLorentzVector lead_quark = bq1.Pt() > bq2.Pt() ? bq1 : bq2;
            double r = lead_quark.Pt()/lead_jet.Pt();
            pdf->Fill(r);
        }
    }

    pdf->Scale(1.0/pdf->Integral());

    auto output = std::make_unique<TFile>("pdf.root","RECREATE");
    pdf->Write();
	output->Write();
	output->Close();

    return 0;
}