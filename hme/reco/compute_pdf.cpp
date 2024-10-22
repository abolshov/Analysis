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

#include "HistManager.hpp"

inline constexpr size_t N_RECO_JETS = 20;
inline constexpr size_t N_BINS = 2000;

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
    auto pdf_mbb = std::make_unique<TH1F>("pdf_mbb", "1d PDF of H->bb mass with true corrections applied", 5*N_BINS, 0, 200);
    auto pdf_numet = std::make_unique<TH1F>("pdf_numet", "1d PDF nu to met pt ratio with true corrections applied", N_BINS, 0, 8);

    HistManager hm;

    std::string nu_met_b("nu_met_b");
    hm.Add(nu_met_b, "nu_met_b", {"ratio", "Count"}, {0, 8}, 100);

    std::string nu_met_b_smear("nu_met_b_smear");
    hm.Add(nu_met_b_smear, "nu_met_b_smear", {"ratio", "Count"}, {0, 8}, 100);

    std::string nu_met_b_smear_lj("nu_met_b_smear_lj");
    hm.Add(nu_met_b_smear_lj, "nu_met_b_smear_lj", {"ratio", "Count"}, {0, 8}, 100);

    std::string nu_met("nu_met");
    hm.Add(nu_met, "nu_met", {"ratio", "Count"}, {0, 8}, 100);

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

        TLorentzVector lead_quark = genb1_p4.Pt() > genb2_p4.Pt() ? genb1_p4 : genb2_p4;
        TLorentzVector subl_quark = genb1_p4.Pt() > genb2_p4.Pt() ? genb2_p4 : genb1_p4;
        double c1 = lead_quark.Pt()/reco_bj1_p4.Pt();
        double c2 = subl_quark.Pt()/reco_bj2_p4.Pt();

        double dpx = -(c1 - 1)*reco_bj1_p4.Px() - (c1 - 1)*reco_bj2_p4.Px();
        double dpy = -(c1 - 1)*reco_bj1_p4.Py() - (c1 - 1)*reco_bj2_p4.Py();

        double dpx_smear = 0.0;
        double dpy_smear = 0.0;
        auto smear_x = std::make_unique<TH1F>("smear_x", "smear_x", 100, -50.0, 50.0);
        auto smear_y = std::make_unique<TH1F>("smear_y", "smear_y", 100, -50.0, 50.0);

        double dpx_lj = 0.0;
        double dpy_lj = 0.0;
        std::vector<std::unique_ptr<TH1F>> smear_lj_px;
        std::vector<std::unique_ptr<TH1F>> smear_lj_py;
        for (size_t i = 0; i < light_jets.size(); ++i)
        {
            TString name_x = Form("smear_lj_px_%d", static_cast<int>(i));
            smear_lj_px.push_back(std::make_unique<TH1F>(name_x, name_x, 100, -50.0, 50.0));

            TString name_y = Form("smear_lj_py_%d", static_cast<int>(i));
            smear_lj_py.push_back(std::make_unique<TH1F>(name_y, name_y, 100, -50.0, 50.0));
        }

        for (int i = 0; i < 1000; ++i)
        {
            smear_x->Fill(rg.Gaus(0, 25.2));
            smear_y->Fill(rg.Gaus(0, 25.2));
            for (size_t i = 0; i < light_jets.size(); ++i)
            {
                auto [dpx, pdy] = GenMETCorrFromJet(rg, light_jets[i], resolutions[i]);
                smear_lj_px[i]->Fill(dpx);
                smear_lj_py[i]->Fill(dpy);
            }
        }

        dpx_smear = smear_x->GetXaxis()->GetBinCenter(smear_x->GetMaximumBin());
        dpy_smear = smear_y->GetXaxis()->GetBinCenter(smear_y->GetMaximumBin());

        for (size_t i = 0; i < light_jets.size(); ++i)
        {
            dpx_lj += smear_lj_px[i]->GetXaxis()->GetBinCenter(smear_lj_px[i]->GetMaximumBin());
            dpy_lj += smear_lj_py[i]->GetXaxis()->GetBinCenter(smear_lj_py[i]->GetMaximumBin());
        }

        {
            TLorentzVector reco_met_corr(reco_met.Px() + dpx, reco_met.Py() + dpy, 0.0, 0.0);
            hm.Fill(nu_met_b, nu.Pt()/reco_met_corr.Pt());
        }

        {
            TLorentzVector reco_met_corr(reco_met.Px() + dpx + dpx_smear, reco_met.Py() + dpy + dpy_smear, 0.0, 0.0);
            hm.Fill(nu_met_b_smear, nu.Pt()/reco_met_corr.Pt());
        }

        {
            TLorentzVector reco_met_corr(reco_met.Px() + dpx + dpx_smear + dpx_lj, reco_met.Py() + dpy + dpy_smear + dpy_lj, 0.0, 0.0);
            hm.Fill(nu_met_b_smear_lj, nu.Pt()/reco_met_corr.Pt());
        }

        TLorentzVector reco_met_corr(reco_met.Px() + dpx + dpx_smear, reco_met.Py() + dpy + dpy_smear, 0.0, 0.0);

        pdf_numet->Fill(nu.Pt()/reco_met_corr.Pt());
        pdf_b1->Fill(c1);
        pdf_b2->Fill(c1);
        pdf_b1b2->Fill(c1, c2); 
        pdf_mbb->Fill((c1*reco_bj1_p4 + c2*reco_bj2_p4).M());

        hm.Fill(nu_met, nu.Pt()/reco_met.Pt());

    }

    hm.DrawStack({nu_met, nu_met_b, nu_met_b_smear, nu_met_b_smear_lj}, "MET PDF comparison", "met_pdf_cmp");

    pdf_b1->Scale(1.0/GetPDFScaleFactor(pdf_b1));
    pdf_b2->Scale(1.0/GetPDFScaleFactor(pdf_b2));
    pdf_mbb->Scale(1.0/GetPDFScaleFactor(pdf_mbb));
    pdf_numet->Scale(1.0/GetPDFScaleFactor(pdf_numet));

    auto output = std::make_unique<TFile>("pdf.root","RECREATE");
    pdf_b1->Write();
    pdf_b2->Write();
    pdf_mbb->Write();
    pdf_numet->Write();
    pdf_b1b2->Write();
	output->Write();
	output->Close();

    return 0;
}