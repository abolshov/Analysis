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

TLorentzVector GenerateResCorrection(TLorentzVector const& v, TRandom3& rg, double res)
{
    TLorentzVector result;
    double dpt = rg.Gaus(0, res);
    double pt = v.Pt();
    while (pt + dpt < 0.0)
    {
        dpt = rg.Gaus(0, res);
    }
    result.SetPtEtaPhiM(pt + dpt, v.Eta(), v.Phi(), v.M());
    return result;
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
    Double_t        genV2_mass;
    Double_t        genV1_mass;

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

    chain->SetBranchAddress("genV2_mass", &genV2_mass);
    chain->SetBranchAddress("genV1_mass", &genV1_mass);

    auto pdf_b1b2 = std::make_unique<TH2F>("pdf_b1b2", "2d PDF simultaneous b jet corrrections", N_BINS, 0, 8, N_BINS, 0, 8);
    auto pdf_hh_dEtadPhi = std::make_unique<TH2F>("pdf_hh_dEtadPhi", "2d PDF dEta vs dPhi between H->bb and H->WW", N_BINS, -8, 8, N_BINS, -8, 8);
    auto pdf_hh_pt_e = std::make_unique<TH2F>("pdf_hh_pt_e", "2d PDF of ratio pt to E of H->bb and H->WW", N_BINS, 0, 1, N_BINS, 0, 1);

    auto pdf_b1 = std::make_unique<TH1F>("pdf_b1", "1d PDF for leading b jet correction", N_BINS, 0, 8);
    auto pdf_b2 = std::make_unique<TH1F>("pdf_b2", "1d PDF for subleading b jet correction", N_BINS, 0, 8);
    auto pdf_mbb = std::make_unique<TH1F>("pdf_mbb", "1d PDF of H->bb mass with true corrections applied", N_BINS, 0, 200);
    auto pdf_numet_pt = std::make_unique<TH1F>("pdf_numet_pt", "1d PDF nu to met pt ratio with true corrections applied", N_BINS, 0, 8);
    auto pdf_numet_dphi = std::make_unique<TH1F>("pdf_numet_dphi", "1d PDF of dPhi between true nu and MET with true corrections applied", N_BINS, -4, 4);
    auto pdf_nulep_deta = std::make_unique<TH1F>("pdf_nulep_deta", "1d PDF of dEta between true nu and reco lep", N_BINS, -8, 8);
    auto pdf_hh_dphi = std::make_unique<TH1F>("pdf_hh_dphi", "1d PDF of dPhi between H->bb and H->WW", N_BINS, -4, 4);
    auto pdf_hh_deta = std::make_unique<TH1F>("pdf_hh_deta", "1d PDF of dEta between H->bb and H->WW", N_BINS, -8, 8);
    auto pdf_mww = std::make_unique<TH1F>("pdf_mww", "1d PDF of H->WW mass with best corrections applied", N_BINS, 0, 200);
    auto pdf_mjj_off = std::make_unique<TH1F>("pdf_mjj_off", "1d PDF of invariant mass of light jets when gen W->qq is offshell", N_BINS, 0, 200);
    auto pdf_mjj_on = std::make_unique<TH1F>("pdf_mjj_on", "1d PDF of invariant mass of light jets when gen W->qq is onshell", N_BINS, 0, 200);
    // auto pdf_mw_had = std::make_unique<TH1F>("pdf_mw_had", "1d PDF of hadronic W mass with best corrections applied", N_BINS, 0, 200);
    // auto pdf_mw_lep = std::make_unique<TH1F>("pdf_mw_lep", "1d PDF of leptonic W mass with reco lep and true nu", N_BINS, 0, 200);
    auto pdf_hbb_pt_e = std::make_unique<TH1F>("pdf_hbb_pt_e", "1d PDF of pt to E ratio for H->bb", N_BINS, 0, 1);
    auto pdf_hww_pt_e = std::make_unique<TH1F>("pdf_hww_pt_e", "1d PDF of pt to E ratio for H->WW", N_BINS, 0, 1);

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
        std::vector<double> resolutions;
        for (int j = 0; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);
            jets.push_back(jet);

            jet *= centralJet_PNetRegPtRawCorr[j];
            resolutions.push_back(jet.Pt()*centralJet_PNetRegPtRawRes[j]);
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

        std::unordered_set<int> match_idx;

        int b1_match = FindBestMatch(genb1_p4, jets);
        int b2_match = FindBestMatch(genb2_p4, jets);
        match_idx.insert(b1_match);
        match_idx.insert(b2_match);

        int q1_match =  FindBestMatch(genq1_p4, jets);
        int q2_match =  FindBestMatch(genq2_p4, jets);
        match_idx.insert(q1_match);
        match_idx.insert(q2_match);

        if (match_idx.size() != 4 || match_idx.count(-1))
        {
            continue;
        }

        TLorentzVector reco_lep_p4;
        reco_lep_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);

        TLorentzVector const& reco_bj1_p4 = jets[b1_match];
        TLorentzVector const& reco_bj2_p4 = jets[b2_match];
        TLorentzVector const& reco_lj1_p4 = jets[q1_match];
        TLorentzVector const& reco_lj2_p4 = jets[q2_match];

        TLorentzVector Hbb_p4, Hww_p4;
        Hbb_p4.SetPtEtaPhiM(genHbb_pt, genHbb_eta, genHbb_phi, genHbb_mass);
        Hww_p4.SetPtEtaPhiM(genHVV_pt, genHVV_eta, genHVV_phi, genHVV_mass);
        pdf_hh_dphi->Fill(Hbb_p4.DeltaPhi(Hww_p4));
        pdf_hh_dEtadPhi->Fill(Hbb_p4.Eta() - Hww_p4.Eta(), Hbb_p4.DeltaPhi(Hww_p4));
        pdf_hh_deta->Fill(Hbb_p4.Eta() - Hww_p4.Eta());

        pdf_hbb_pt_e->Fill(Hbb_p4.Pt()/Hbb_p4.E());
        pdf_hww_pt_e->Fill(Hww_p4.Pt()/Hww_p4.E());
        pdf_hh_pt_e->Fill(Hbb_p4.Pt()/Hbb_p4.E(), Hww_p4.Pt()/Hww_p4.E());

        TLorentzVector reco_met;
        reco_met.SetPtEtaPhiM(PuppiMET_pt, 0.0, PuppiMET_phi, 0.0);

        TLorentzVector nu;
        nu.SetPtEtaPhiM(genV1prod2_pt, genV1prod2_eta, genV1prod2_phi, genV1prod2_mass);

        double c1 = genb1_p4.Pt()/reco_bj1_p4.Pt();
        double c2 = genb2_p4.Pt()/reco_bj2_p4.Pt();

        double dpx_b = -(c1 - 1)*reco_bj1_p4.Px() - (c2 - 1)*reco_bj2_p4.Px();
        double dpy_b = -(c1 - 1)*reco_bj1_p4.Py() - (c2 - 1)*reco_bj2_p4.Py();

        double dpx_smear = 0.0;
        double dpy_smear = 0.0;
        auto smear_x = std::make_unique<TH1F>("smear_x", "smear_x", 100, -100.0, 100.0);
        auto smear_y = std::make_unique<TH1F>("smear_y", "smear_y", 100, -100.0, 100.0);

        double dpx_l = 0.0;
        double dpy_l = 0.0;
        auto res_light_x = std::make_unique<TH1F>("res_light_x", "res_light_x", 100, -100.0, 100.0);
        auto res_light_y = std::make_unique<TH1F>("res_light_y", "res_light_y", 100, -100.0, 100.0);

        auto hww_mass_hist = std::make_unique<TH1F>("hww_mass_hist", "hww_mass_hist", 100, 0, 300.0);

        double mass_diff = 10e4;
        TLorentzVector reco_Wlep_p4 = reco_lep_p4 + nu;
        TLorentzVector reco_Whad_p4;
        for (int i = 0; i < 1000; ++i)
        {
            smear_x->Fill(rg.Gaus(0, 25.2));
            smear_y->Fill(rg.Gaus(0, 25.2));

            TLorentzVector lj1_corr_p4 = GenerateResCorrection(reco_lj1_p4, rg, resolutions[q1_match]);
            TLorentzVector lj2_corr_p4 = GenerateResCorrection(reco_lj2_p4, rg, resolutions[q2_match]);

            hww_mass_hist->Fill((reco_Wlep_p4 + lj1_corr_p4 + lj2_corr_p4).M());

            double dm = std::abs((lj1_corr_p4 + lj2_corr_p4).M() - genV2_mass);
            if (dm < mass_diff)
            {
                reco_Whad_p4 = lj1_corr_p4 + lj2_corr_p4;
            }
            mass_diff = std::min(dm, mass_diff);

            res_light_x->Fill(lj1_corr_p4.Px() + lj2_corr_p4.Px() - reco_lj1_p4.Px() - reco_lj2_p4.Px());
            res_light_x->Fill(lj1_corr_p4.Py() + lj2_corr_p4.Py() - reco_lj1_p4.Py() - reco_lj2_p4.Py());
        }
        TLorentzVector reco_Hww_p4 = reco_Whad_p4 + reco_Wlep_p4;
        // pdf_mw_had->Fill(reco_Whad_p4.M());
        pdf_mww->Fill(reco_Hww_p4.M());

        // double hww_mass = hww_mass_hist->GetXaxis()->GetBinCenter(hww_mass_hist->GetMaximumBin());
        // double hww_mass = (reco_Wlep_p4 + reco_lj1_p4 + reco_lj2_p4).M();
        // pdf_mww->Fill(hww_mass);

        if (genV1_mass > genV2_mass)
        {
            pdf_mjj_off->Fill((reco_lj1_p4 + reco_lj2_p4).M());
        }
        else 
        {
            pdf_mjj_on->Fill((reco_lj1_p4 + reco_lj2_p4).M());
        }

        dpx_smear = smear_x->GetXaxis()->GetBinCenter(smear_x->GetMaximumBin());
        dpy_smear = smear_y->GetXaxis()->GetBinCenter(smear_y->GetMaximumBin());

        dpx_l = res_light_x->GetXaxis()->GetBinCenter(res_light_x->GetMaximumBin());
        dpy_l = res_light_x->GetXaxis()->GetBinCenter(res_light_x->GetMaximumBin());

        double met_corr_px = reco_met.Px() + dpx_b + dpx_smear - dpx_l;
        double met_corr_py = reco_met.Py() + dpy_b + dpy_smear - dpy_l;
        double met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
        double met_corr_phi = std::atan2(met_corr_py, met_corr_px);
        TLorentzVector reco_met_corr;
        reco_met_corr.SetPtEtaPhiM(met_corr_pt, 0.0, met_corr_phi, 0.0);

        pdf_numet_pt->Fill(nu.Pt()/reco_met_corr.Pt());
        pdf_numet_dphi->Fill(nu.DeltaPhi(reco_met_corr));
        pdf_b1->Fill(c1);
        pdf_b2->Fill(c2);
        pdf_b1b2->Fill(c1, c2); 
        pdf_mbb->Fill((c1*reco_bj1_p4 + c2*reco_bj2_p4).M());
        pdf_nulep_deta->Fill(nu.Eta() - reco_lep_p4.Eta());
    }

    pdf_b1->Scale(1.0/GetPDFScaleFactor(pdf_b1));
    pdf_b2->Scale(1.0/GetPDFScaleFactor(pdf_b2));
    pdf_mbb->Scale(1.0/GetPDFScaleFactor(pdf_mbb));
    pdf_numet_pt->Scale(1.0/GetPDFScaleFactor(pdf_numet_pt));
    pdf_numet_dphi->Scale(1.0/GetPDFScaleFactor(pdf_numet_dphi));
    pdf_hh_dphi->Scale(1.0/GetPDFScaleFactor(pdf_hh_dphi));
    pdf_hh_deta->Scale(1.0/GetPDFScaleFactor(pdf_hh_deta));
    pdf_nulep_deta->Scale(1.0/GetPDFScaleFactor(pdf_nulep_deta));
    pdf_mww->Scale(1.0/GetPDFScaleFactor(pdf_mww));
    // pdf_mw_had->Scale(1.0/GetPDFScaleFactor(pdf_mw_had));
    pdf_mjj_off->Scale(1.0/GetPDFScaleFactor(pdf_mjj_off));
    pdf_mjj_on->Scale(1.0/GetPDFScaleFactor(pdf_mjj_on));
    pdf_b1b2->Scale(1.0/GetPDFScaleFactor(pdf_b1b2));
    pdf_hh_dEtadPhi->Scale(1.0/GetPDFScaleFactor(pdf_hh_dEtadPhi));
    pdf_hh_pt_e->Scale(1.0/GetPDFScaleFactor(pdf_hh_pt_e));
    pdf_hbb_pt_e->Scale(1.0/GetPDFScaleFactor(pdf_hbb_pt_e));
    pdf_hww_pt_e->Scale(1.0/GetPDFScaleFactor(pdf_hww_pt_e));

    auto output = std::make_unique<TFile>("pdf.root", "RECREATE");
    pdf_b1->Write();
    pdf_b2->Write();
    pdf_mbb->Write();
    pdf_hh_dphi->Write();
    pdf_numet_pt->Write();
    pdf_numet_dphi->Write();
    pdf_hh_deta->Write();
    pdf_b1b2->Write();
    pdf_nulep_deta->Write();
    pdf_mww->Write();
    pdf_hh_dEtadPhi->Write();
    pdf_mjj_off->Write();
    pdf_mjj_on->Write();
    pdf_hh_pt_e->Write();
    pdf_hbb_pt_e->Write();
    pdf_hww_pt_e->Write();
    // pdf_mw_had->Write();
	output->Write();
	output->Close();

    return 0;
}