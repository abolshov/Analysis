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
#include "TRandom3.h"

inline constexpr size_t N_RECO_JETS = 20;
inline constexpr size_t N_BINS = 100;

template <class T>
double GetPDFScaleFactor(std::unique_ptr<T> const& hist)
{
    int binmax = hist->GetMaximumBin();
    return hist->GetBinContent(binmax);
}

// generates dpx, dpy due to jet PNet resolution
std::pair<double, double> GenMETCorrFromJet(TRandom3& rg, TLorentzVector p, double res)
{
    double dpt = rg.Gaus(0, res);
    bool valid = p.Pt() + dpt > 0.0;
    double dpx = valid ? p.Px() : 0.0;
    double dpy = valid ? p.Py() : 0.0;
    if (valid)
    {
        p.SetPtEtaPhiM(p.Pt() + dpt, p.Eta(), p.Phi(), p.M());
        
        dpx -= p.Px();
        dpx *= -1.0;
        
        dpy -= p.Py();
        dpy *= -1.0;
    }
    return {dpx, dpy};   
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

    chain->SetBranchAddress("lep1_pt", &lep1_pt);
    chain->SetBranchAddress("lep1_eta", &lep1_eta);
    chain->SetBranchAddress("lep1_phi", &lep1_phi);
    chain->SetBranchAddress("lep1_mass", &lep1_mass);

    chain->SetBranchAddress("lep1_type", &lep1_type);
    chain->SetBranchAddress("lep1_genLep_kind", &lep1_genLep_kind);

    Double_t        genHVV_pt;
    Double_t        genHVV_eta;
    Double_t        genHVV_phi;
    Double_t        genHVV_mass;    

    Double_t        genHbb_pt;
    Double_t        genHbb_eta;
    Double_t        genHbb_phi;
    Double_t        genHbb_mass;

    chain->SetBranchAddress("genHVV_pt", &genHVV_pt);
    chain->SetBranchAddress("genHVV_eta", &genHVV_eta);
    chain->SetBranchAddress("genHVV_phi", &genHVV_phi);
    chain->SetBranchAddress("genHVV_mass", &genHVV_mass);

    chain->SetBranchAddress("genHbb_pt", &genHbb_pt);
    chain->SetBranchAddress("genHbb_eta", &genHbb_eta);
    chain->SetBranchAddress("genHbb_phi", &genHbb_phi);
    chain->SetBranchAddress("genHbb_mass", &genHbb_mass);    

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

    Float_t         PuppiMET_pt;
    Float_t         PuppiMET_phi;

    // neutrino
    Double_t        genV1prod2_pt;
    Double_t        genV1prod2_eta;
    Double_t        genV1prod2_phi;
    Double_t        genV1prod2_mass;

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

    auto pdf_b1 = std::make_unique<TH1F>("pdf_b1", "1d PDF for leading b jet correction", N_BINS, 0, 8);
    auto pdf_b2 = std::make_unique<TH1F>("pdf_b2", "1d PDF for subleading b jet correction", N_BINS, 0, 8);
    auto pdf_b1b2 = std::make_unique<TH2F>("pdf_b1b2", "2d PDF simultaneous b jet corrrections", N_BINS, 0, 8, N_BINS, 0, 8);
    auto pdf_mbb = std::make_unique<TH1F>("pdf_mbb", "1d PDF of H->bb mass with true corrections applied", N_BINS, 0, 200);
    auto pdf_numet_pt = std::make_unique<TH1F>("pdf_numet_pt", "1d PDF nu to met pt ratio with true corrections applied", N_BINS, 0, 8);
    auto pdf_numet_dphi = std::make_unique<TH1F>("pdf_numet_dphi", "1d PDF of dPhi between true nu and MET with true corrections applied", N_BINS, -4, 4);
    auto pdf_nulep_deta = std::make_unique<TH1F>("pdf_nulep_deta", "1d PDF of dEta between true nu and reco lep", N_BINS, -8, 8);
    auto pdf_hh_dphi = std::make_unique<TH1F>("pdf_hh_dphi", "1d PDF of dPhi between H->bb and H->WW", N_BINS, -4, 4);

    TRandom3 rg;
    rg.SetSeed(42);

    long long nEvents = chain->GetEntries();
    for (long long i = 0; i < nEvents; ++i)
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

        std::vector<TLorentzVector> jets;
        for (int j = 0; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);
            jets.push_back(jet);
        }

        int b1_match = FindBestMatch(genb1_p4, jets);
        int b2_match = FindBestMatch(genb2_p4, jets);
        if (b1_match == -1 || b2_match == -1 || b1_match == b2_match)
        {
            continue;
        }

        bool reco_lep1_mu = (lep1_type == 2);
        bool reco_lep1_ele = (lep1_type == 1);

        bool gen_lep1_mu = ((lep1_genLep_kind == 2) || (lep1_genLep_kind == 4));
        bool gen_lep1_ele = ((lep1_genLep_kind == 1) || (lep1_genLep_kind == 3));

        bool corr_lep_reco = ((reco_lep1_mu && gen_lep1_mu) || (reco_lep1_ele && gen_lep1_ele));
        
        if (!corr_lep_reco)
        {
            continue;
        }

        TLorentzVector reco_bj1_p4 = jets[b1_match];
        TLorentzVector reco_bj2_p4 = jets[b2_match];

        TLorentzVector Hbb_p4, Hww_p4;
        Hbb_p4.SetPtEtaPhiM(genHbb_pt, genHbb_eta, genHbb_phi, genHbb_mass);
        Hww_p4.SetPtEtaPhiM(genHVV_pt, genHVV_eta, genHVV_phi, genHVV_mass);
        pdf_hh_dphi->Fill(Hbb_p4.DeltaPhi(Hww_p4));

        std::vector<TLorentzVector> light_jets;
        std::vector<double> resolutions;
        for (int j = 2; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);

            // use PNet correction to correct p4 of light jets
            jet *= centralJet_PNetRegPtRawCorr[j];
            light_jets.push_back(jet);

            // save resolutions of pt corrections of light jets
            resolutions.push_back(centralJet_pt[j]*centralJet_PNetRegPtRawRes[j]);
        }

        TLorentzVector reco_met;
        reco_met.SetPtEtaPhiM(PuppiMET_pt, 0.0, PuppiMET_phi, 0.0);

        TLorentzVector nu;
        nu.SetPtEtaPhiM(genV1prod2_pt, genV1prod2_eta, genV1prod2_phi, genV1prod2_mass);

        double c1 = genb1_p4.Pt()/reco_bj1_p4.Pt();
        double c2 = genb2_p4.Pt()/reco_bj2_p4.Pt();

        double dpx = -(c1 - 1)*reco_bj1_p4.Px() - (c1 - 1)*reco_bj2_p4.Px();
        double dpy = -(c1 - 1)*reco_bj1_p4.Py() - (c1 - 1)*reco_bj2_p4.Py();

        double dpx_smear = 0.0;
        double dpy_smear = 0.0;
        auto smear_x = std::make_unique<TH1F>("smear_x", "smear_x", 100, -50.0, 50.0);
        auto smear_y = std::make_unique<TH1F>("smear_y", "smear_y", 100, -50.0, 50.0);

        for (int i = 0; i < 1000; ++i)
        {
            smear_x->Fill(rg.Gaus(0, 25.2));
            smear_y->Fill(rg.Gaus(0, 25.2));
        }

        dpx_smear = smear_x->GetXaxis()->GetBinCenter(smear_x->GetMaximumBin());
        dpy_smear = smear_y->GetXaxis()->GetBinCenter(smear_y->GetMaximumBin());

        TLorentzVector reco_met_corr(reco_met.Px() + dpx + dpx_smear, reco_met.Py() + dpy + dpy_smear, 0.0, 0.0);

        pdf_numet_pt->Fill(nu.Pt()/reco_met_corr.Pt());
        pdf_numet_dphi->Fill(nu.DeltaPhi(reco_met_corr));
        pdf_b1->Fill(c1);
        pdf_b2->Fill(c2);
        pdf_b1b2->Fill(c1, c2); 
        pdf_mbb->Fill((c1*reco_bj1_p4 + c2*reco_bj2_p4).M());

        TLorentzVector reco_lep_p4;
        reco_lep_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
        pdf_nulep_deta->Fill(nu.Eta() - reco_lep_p4.Eta());
    }

    pdf_b1->Scale(1.0/GetPDFScaleFactor(pdf_b1));
    pdf_b2->Scale(1.0/GetPDFScaleFactor(pdf_b2));
    pdf_mbb->Scale(1.0/GetPDFScaleFactor(pdf_mbb));
    pdf_numet_pt->Scale(1.0/GetPDFScaleFactor(pdf_numet_pt));
    pdf_numet_dphi->Scale(1.0/GetPDFScaleFactor(pdf_numet_dphi));
    pdf_hh_dphi->Scale(1.0/GetPDFScaleFactor(pdf_hh_dphi));
    pdf_nulep_deta->Scale(1.0/GetPDFScaleFactor(pdf_nulep_deta));

    auto output = std::make_unique<TFile>("pdf.root", "RECREATE");
    pdf_b1->Write();
    pdf_b2->Write();
    pdf_mbb->Write();
    pdf_hh_dphi->Write();
    pdf_numet_pt->Write();
    pdf_numet_dphi->Write();
    pdf_b1b2->Write();
    pdf_nulep_deta->Write();
	output->Write();
	output->Close();

    return 0;
}