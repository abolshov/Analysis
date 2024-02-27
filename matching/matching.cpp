#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include "MatchingTools.hpp"
#include "HistManager.hpp"

// static constexpr int MAX_GENJET = 16;
static constexpr int MAX_GENJET = 20;
static constexpr int MAX_GENPART = 270;

void save_2d_dist(TH2F* dist, 
                  std::string const& name,
                  std::string const& title_x,
                  std::string const& title_y)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();
    dist->SetStats(0);
    TStyle* gStyle = new TStyle();
    gStyle->SetPalette(kRainBow);

    dist->GetXaxis()->SetTitle(title_x.c_str());
    dist->GetYaxis()->SetTitle(title_y.c_str());
    dist->SetLineWidth(3);
    dist->Draw("colz");
    c1->SaveAs((name + ".png").c_str());

    delete gStyle;
    delete c1;
}

static constexpr double MATCH_PT_THRESH = 3.0;

int main()
{
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    // TFile *myFile = TFile::Open("NanoAOD_600GeV_1000Events.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    UInt_t          nGenPart;
    Float_t         GenPart_eta[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_mass[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_phi[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_pt[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_pdgId[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_status[MAX_GENPART];   //[nGenPart]
    // Int_t           GenPart_statusFlags[MAX_GENPART];   //[nGenPart]

    UInt_t          nGenJet;
    Float_t         GenJet_eta[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_mass[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_phi[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_pt[MAX_GENJET];   //[nGenJet]

    Int_t           GenJet_partonFlavour[MAX_GENJET];   //[nGenJet]
    // UChar_t         GenJet_hadronFlavour[MAX_GENJET];   //[nGenJet] doesn't work for some reason - prints empty spaces
    Int_t           Jet_genJetIdx[MAX_GENJET];  

    myTree->SetBranchAddress("nGenPart", &nGenPart);
    myTree->SetBranchAddress("GenPart_eta", &GenPart_eta);
    myTree->SetBranchAddress("GenPart_mass", &GenPart_mass);
    myTree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    myTree->SetBranchAddress("GenPart_pt", &GenPart_pt);
    myTree->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    myTree->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    myTree->SetBranchAddress("GenPart_status", &GenPart_status);
    myTree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    // myTree->SetBranchAddress("GenPart_statusFlags", &GenPart_statusFlags);

    myTree->SetBranchAddress("nGenJet", &nGenJet);
    myTree->SetBranchAddress("GenJet_eta", &GenJet_eta);
    myTree->SetBranchAddress("GenJet_mass", &GenJet_mass);
    myTree->SetBranchAddress("GenJet_phi", &GenJet_phi);
    myTree->SetBranchAddress("GenJet_pt", &GenJet_pt);

    myTree->SetBranchAddress("GenJet_partonFlavour", &GenJet_partonFlavour);
    myTree->SetBranchAddress("Jet_genJetIdx", &Jet_genJetIdx);

    int nEvents = myTree->GetEntries();

    // counters
    int non_empty_sig = 0;
    int valid_sig = 0;
    int matching_fails = 0;
    int accep_evts = 0;

    int n_q1_exc = 0;
    int n_q2_exc = 0;
    int n_any_exc = 0;

    // hists
    HistManager hm(10);

    int nbins = 101;
    std::pair<double, double> pt_ratio_range{0.0, 6.0};
    std::pair<std::string, std::string> mass_labels{"[GeV]", "AU"};
    std::pair<std::string, std::string> pt_ratio_label_labels{"quark pt/jet pt", "AU"};

    hm.Add("w_from_quarks", "W mass from quarks", mass_labels, {0, 120}, nbins);
    hm.Add("w_from_jets", "W mass from jets", mass_labels, {0, 150}, nbins/2);
    hm.Add("h_from_quarks", "Higgs mass from quarks", mass_labels, {120, 130}, nbins);
    hm.Add("h_from_jets", "Higgs mass from jets", mass_labels, {50, 175}, nbins/2);
    hm.Add("lead_b_pt_ratio_hist", "leading b pair pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);
    hm.Add("lead_l_pt_ratio_hist", "leading light pair pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);
    hm.Add("sublead_b_pt_ratio_hist", "subleading b pair pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);
    hm.Add("sublead_l_pt_ratio_hist", "subleading light pair pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);
    hm.Add("overlap_pt_ratio_1", "pair 1 pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);
    hm.Add("overlap_pt_ratio_2", "pair 2 pt ratio", pt_ratio_label_labels, pt_ratio_range, nbins);

    std::string w_from_quarks("w_from_quarks");
    std::string w_from_jets("w_from_jets");
    std::string h_from_quarks("h_from_quarks");
    std::string h_from_jets("h_from_jets");
    std::string lead_b_pt_ratio_hist("lead_b_pt_ratio_hist");
    std::string lead_l_pt_ratio_hist("lead_l_pt_ratio_hist");
    std::string sublead_b_pt_ratio_hist("sublead_b_pt_ratio_hist");
    std::string sublead_l_pt_ratio_hist("sublead_l_pt_ratio_hist");
    std::string overlap_pt_ratio_1("overlap_pt_ratio_1");
    std::string overlap_pt_ratio_2("overlap_pt_ratio_2");

    std::cout << std::boolalpha;

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);
        if (!sig.empty())
        {
            ++non_empty_sig;
            if (CheckSignal(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                ++valid_sig;

                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet{GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, static_cast<int>(nGenJet)};

                int b_idx = sig[SIG::b];
                int bbar_idx = sig[SIG::bbar];
                int q1_idx = sig[SIG::q1];
                int q2_idx = sig[SIG::q2];
                
                int b_match = Match(b_idx, genpart, genjet);
                int bbar_match = Match(bbar_idx, genpart, genjet);
                int q1_match = Match(q1_idx, genpart, genjet);
                int q2_match = Match(q2_idx, genpart, genjet);

                // skip event if matching failed
                std::vector<int> matches{b_match, bbar_match, q1_match, q2_match};
                if(std::any_of(matches.begin(), matches.end(), [](int x) { return x == -1;}) || std::unique(matches.begin(), matches.end()) != matches.end())
                {
                    ++matching_fails;
                    continue;
                } 
                else
                {
                    ++accep_evts;

                    TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                    TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);
                    TLorentzVector bj_p4 = GetP4(genjet, b_match);
                    TLorentzVector bbarj_p4 = GetP4(genjet, bbar_match);

                    TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                    TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                    TLorentzVector lj1_p4 = GetP4(genjet, q1_match);
                    TLorentzVector lj2_p4 = GetP4(genjet, q2_match);

                    if (lj1_p4.DeltaR(lj2_p4) < 2*DR_THRESH)
                    {
                        hm.Fill(overlap_pt_ratio_1, lq1_p4.Pt()/lj1_p4.Pt());
                        hm.Fill(overlap_pt_ratio_2, lq2_p4.Pt()/lj2_p4.Pt());
                    }

                    TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                    double dR_lj1_lep = lj1_p4.DeltaR(l_p4);
                    double dR_lj2_lep = lj2_p4.DeltaR(l_p4);
                    if (dR_lj1_lep < DR_THRESH && dR_lj2_lep > DR_THRESH)
                    {
                        TLorentzVector tmp = lj1_p4 - l_p4;
                        double resc = tmp.Pt()/lj1_p4.Pt();
                        lj1_p4 *= resc;
                    }
                    if (dR_lj2_lep < DR_THRESH && dR_lj1_lep > DR_THRESH) 
                    {
                        TLorentzVector tmp = lj2_p4 - l_p4;
                        double resc = tmp.Pt()/lj2_p4.Pt();
                        lj2_p4 *= resc;
                    }

                    double leading_b_pt_ratio = (bq_p4.Pt() > bbarq_p4.Pt()) ? (bq_p4.Pt()/bj_p4.Pt()) : (bbarq_p4.Pt()/bbarj_p4.Pt());
                    double leading_l_pt_ratio = (lq1_p4.Pt() > lq2_p4.Pt()) ? (lq1_p4.Pt()/lj1_p4.Pt()) : (lq2_p4.Pt()/lj2_p4.Pt());

                    double subleading_b_pt_ratio = (bq_p4.Pt() > bbarq_p4.Pt()) ? (bbarq_p4.Pt()/bbarj_p4.Pt()) : (bq_p4.Pt()/bj_p4.Pt());
                    double subleading_l_pt_ratio = (lq1_p4.Pt() > lq2_p4.Pt()) ? (lq2_p4.Pt()/lj2_p4.Pt()) : (lq1_p4.Pt()/lj1_p4.Pt());

                    hm.Fill(lead_b_pt_ratio_hist, leading_b_pt_ratio);
                    hm.Fill(lead_l_pt_ratio_hist, leading_l_pt_ratio);
                    hm.Fill(sublead_b_pt_ratio_hist, subleading_b_pt_ratio);
                    hm.Fill(sublead_l_pt_ratio_hist, subleading_l_pt_ratio);

                    #ifdef DEBUG
                        std::vector<double> ratios{ leading_b_pt_ratio, subleading_b_pt_ratio, leading_l_pt_ratio, subleading_l_pt_ratio };
                        if (std::any_of(ratios.begin(), ratios.end(), [](double x) { return x > MATCH_PT_THRESH; }))
                        {
                            MatchIndex mi = {{b_idx, b_match}, {bbar_idx, bbar_match}, {q1_idx, q1_match}, {q2_idx, q2_match}};
                            MatchKinematics mk = {genpart, genjet};
                            std::pair<int*, int*> ptrs = {GenPart_genPartIdxMother, GenPart_pdgId};
                            DrawEventMap(mk, mi, i, ptrs);
                        }
                    #endif

                    hm.Fill(w_from_quarks, (lq1_p4 + lq2_p4).M());
                    hm.Fill(w_from_jets, (lj1_p4 + lj2_p4).M());
                    hm.Fill(h_from_quarks, (bq_p4 + bbarq_p4).M());
                    hm.Fill(h_from_jets, (bj_p4 + bbarj_p4).M());
                }
            }
        }
    }

    std::cout << "Finished processing\n";
    hm.Draw();
    hm.DrawStack({w_from_quarks, w_from_jets}, "W mass", "w_mass");
    hm.DrawStack({h_from_quarks, h_from_jets}, "Higgs mass", "higgs_mass");
    hm.DrawStack({lead_b_pt_ratio_hist, lead_l_pt_ratio_hist, sublead_b_pt_ratio_hist, sublead_l_pt_ratio_hist}, "Pt ratios", "pt_ratios");
    hm.DrawStack({lead_b_pt_ratio_hist, lead_l_pt_ratio_hist}, "leading pair Pt ratios", "lead_pt_ratios");
    hm.DrawStack({sublead_b_pt_ratio_hist, sublead_l_pt_ratio_hist}, "subleading pair Pt ratios", "sublead_pt_ratios");
    hm.DrawStack({lead_b_pt_ratio_hist, sublead_b_pt_ratio_hist}, "bb pair Pt ratios", "bb_pt_ratios");
    hm.DrawStack({lead_l_pt_ratio_hist, sublead_l_pt_ratio_hist}, "light pair Pt ratios", "qq_pt_ratios");
    hm.DrawStack({overlap_pt_ratio_1, overlap_pt_ratio_2}, "overlapping light pair Pt ratios", "overlap_qq_pt_ratios");

    std::cout << "nEvents = " << nEvents << "\n" 
              << "\tnon_empty_sig = " << non_empty_sig << "\n" 
              << "\tvalid_sig = " << valid_sig << "\n"
              << "\tmatching_fails = " << matching_fails << "\n"
              << "\taccep_evts = " << accep_evts << "\n"
              << "\tn_q1_exc = " << n_q1_exc << "\n"
              << "\tn_q2_exc = " << n_q2_exc << "\n"
              << "\tn_any_exc = " << n_any_exc << "\n";

    return 0;
}