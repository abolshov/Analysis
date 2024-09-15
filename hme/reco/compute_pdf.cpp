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
#include "TCanvas.h"

static constexpr int N_RECO_JETS = 20;

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

    // V2: hadronic
    Double_t        genV2prod1_pt;
    Double_t        genV2prod1_vis_pt;
    Double_t        genV2prod1_eta;
    Double_t        genV2prod1_vis_eta;
    Double_t        genV2prod1_phi;
    Double_t        genV2prod1_vis_phi;
    Double_t        genV2prod1_mass;
    Double_t        genV2prod1_vis_mass;

    Double_t        genV2prod2_pt;
    Double_t        genV2prod2_vis_pt;
    Double_t        genV2prod2_eta;
    Double_t        genV2prod2_vis_eta;
    Double_t        genV2prod2_phi;
    Double_t        genV2prod2_vis_phi;
    Double_t        genV2prod2_mass;
    Double_t        genV2prod2_vis_mass;

    // H->bb
    Double_t        genb1_pt;
    Double_t        genb1_vis_pt;
    Double_t        genb1_eta;
    Double_t        genb1_vis_eta;
    Double_t        genb1_phi;
    Double_t        genb1_vis_phi;
    Double_t        genb1_mass;
    Double_t        genb1_vis_mass;
    Double_t        genb2_pt;
    Double_t        genb2_vis_pt;
    Double_t        genb2_eta;
    Double_t        genb2_vis_eta;
    Double_t        genb2_phi;
    Double_t        genb2_vis_phi;
    Double_t        genb2_mass;
    Double_t        genb2_vis_mass;

    chain->SetBranchAddress("ncentralJet", &ncentralJet);
    chain->SetBranchAddress("centralJet_pt", centralJet_pt);
    chain->SetBranchAddress("centralJet_eta", centralJet_eta);
    chain->SetBranchAddress("centralJet_phi", centralJet_phi);
    chain->SetBranchAddress("centralJet_mass", centralJet_mass);

    chain->SetBranchAddress("genb1_pt", &genb1_pt);
    chain->SetBranchAddress("genb1_eta", &genb1_eta);
    chain->SetBranchAddress("genb1_phi", &genb1_phi);
    chain->SetBranchAddress("genb1_mass", &genb1_mass);

    chain->SetBranchAddress("genb2_pt", &genb2_pt);
    chain->SetBranchAddress("genb2_eta", &genb2_eta);
    chain->SetBranchAddress("genb2_phi", &genb2_phi);
    chain->SetBranchAddress("genb2_mass", &genb2_mass);

    chain->SetBranchAddress("genV2prod1_pt", &genV2prod1_pt);
    chain->SetBranchAddress("genV2prod1_vis_pt", &genV2prod1_vis_pt);
    chain->SetBranchAddress("genV2prod1_eta", &genV2prod1_eta);
    chain->SetBranchAddress("genV2prod1_vis_eta", &genV2prod1_vis_eta);
    chain->SetBranchAddress("genV2prod1_phi", &genV2prod1_phi);
    chain->SetBranchAddress("genV2prod1_vis_phi", &genV2prod1_vis_phi);
    chain->SetBranchAddress("genV2prod1_mass", &genV2prod1_mass);
    chain->SetBranchAddress("genV2prod1_vis_mass", &genV2prod1_vis_mass);

    chain->SetBranchAddress("genV2prod2_pt", &genV2prod2_pt);
    chain->SetBranchAddress("genV2prod2_vis_pt", &genV2prod2_vis_pt);
    chain->SetBranchAddress("genV2prod2_eta", &genV2prod2_eta);
    chain->SetBranchAddress("genV2prod2_vis_eta", &genV2prod2_vis_eta);
    chain->SetBranchAddress("genV2prod2_phi", &genV2prod2_phi);
    chain->SetBranchAddress("genV2prod2_vis_phi", &genV2prod2_vis_phi);
    chain->SetBranchAddress("genV2prod2_mass", &genV2prod2_mass);
    chain->SetBranchAddress("genV2prod2_vis_mass", &genV2prod2_vis_mass);

    chain->SetBranchAddress("genb1_pt", &genb1_pt);
    chain->SetBranchAddress("genb1_vis_pt", &genb1_vis_pt);
    chain->SetBranchAddress("genb1_eta", &genb1_eta);
    chain->SetBranchAddress("genb1_vis_eta", &genb1_vis_eta);
    chain->SetBranchAddress("genb1_phi", &genb1_phi);
    chain->SetBranchAddress("genb1_vis_phi", &genb1_vis_phi);
    chain->SetBranchAddress("genb1_mass", &genb1_mass);
    chain->SetBranchAddress("genb1_vis_mass", &genb1_vis_mass);
    chain->SetBranchAddress("genb2_pt", &genb2_pt);
    chain->SetBranchAddress("genb2_vis_pt", &genb2_vis_pt);
    chain->SetBranchAddress("genb2_eta", &genb2_eta);
    chain->SetBranchAddress("genb2_vis_eta", &genb2_vis_eta);
    chain->SetBranchAddress("genb2_phi", &genb2_phi);
    chain->SetBranchAddress("genb2_vis_phi", &genb2_vis_phi);
    chain->SetBranchAddress("genb2_mass", &genb2_mass);
    chain->SetBranchAddress("genb2_vis_mass", &genb2_vis_mass);

    auto pdf = std::make_unique<TH1F>("pdf", "pdf", 1000, 0, 6);

    int nEvents = chain->GetEntries();
    for (int i = 0; i < nEvents; ++i)
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

        // filter out resolved events
        if (genb1_p4.DeltaR(genb2_p4) < 0.4 || genq1_p4.DeltaR(genq2_p4) < 0.4)
        {
            continue;
        }

        // require matching for all jets
        bool all_matched = (genb1_vis_pt > 0.0) && (genb2_vis_pt > 0.0) && (genV2prod1_vis_pt > 0.0) && (genV2prod2_vis_pt); 
        if (!all_matched)
        {
            continue;
        }

        // TLorentzVector gen_bj1_p4, gen_bj2_p4, gen_lj1_p4, gen_lj2_p4; // gen jets matched to quarks
        // gen_bj1_p4.SetPtEtaPhiM(genb1_vis_pt, genb1_vis_eta, genb1_vis_phi, genb1_vis_mass);
        // gen_bj2_p4.SetPtEtaPhiM(genb2_vis_pt, genb2_vis_eta, genb2_vis_phi, genb2_vis_mass);
        // gen_lj1_p4.SetPtEtaPhiM(genV2prod1_vis_pt, genV2prod1_vis_eta, genV2prod1_vis_phi, genV2prod1_vis_mass);
        // gen_lj2_p4.SetPtEtaPhiM(genV2prod2_vis_pt, genV2prod2_vis_eta, genV2prod2_vis_phi, genV2prod2_vis_mass);

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
        double r = lead_quark.Pt()/reco_bj1_p4.Pt();
        pdf->Fill(r);
    }

    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    pdf->SetLineWidth(2);
    pdf->Draw("hist");
    pdf->Scale(1.0/pdf->Integral());

    c1->SaveAs("pdf.png");

    auto output = std::make_unique<TFile>("pdf.root","RECREATE");
    pdf->Write();
	output->Write();
	output->Close();

    return 0;
}