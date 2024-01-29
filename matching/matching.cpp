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

// static constexpr int MAX_GENJET = 16;
static constexpr int MAX_GENJET = 20;
static constexpr int MAX_GENPART = 270;

void save_1d_dist(TH1F* dist, std::string const& name, std::string const& title)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

    dist->GetXaxis()->SetTitle(title.c_str());
    dist->SetLineWidth(3);
    dist->SetStats(1);
    dist->DrawNormalized();
    c1->SaveAs((name + ".png").c_str());

    delete c1;
}

void save_1d_stack(std::vector<TH1F*> const& distrs,
                   std::vector<std::string> const& legends,
                   std::string const& name,
                   std::string const& title,
                   std::string const& axis_label)
{
    if (distrs.size() != legends.size())
    {
        std::cout << "number of legends and histograms do not match!" << std::endl;
        return;
    }
    TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

    THStack* stack = new THStack("stack", title.c_str());
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (int i = 0; i < static_cast<int>(distrs.size()); ++i)
    {
        distrs[i]->SetLineWidth(3);
        int line_color = i + 1 < 5 ? i + 1 : i + 2;
        distrs[i]->SetLineColor(line_color);
        stack->Add(distrs[i]);
        legend->AddEntry(distrs[i], legends[i].c_str(), "l");
    }

    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle(axis_label.c_str());
    legend->Draw();
    c1->SaveAs((name + ".png").c_str());

    delete stack;
    delete legend;
    delete c1;
}

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
    TStyle* gStyle = new TStyle();
    gStyle->SetPalette(kRainBow);
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    TFile *myFile = TFile::Open("NanoAOD_600GeV_1000Events.root");
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

    int nEvents = myTree->GetEntries();

    // counters
    int non_empty_sig = 0;
    int valid_sig = 0;
    int matching_fails = 0;
    int fill_cnt = 0;

    // hists
    int nbins = 101;
    auto lead_b_pt_ratio_hist = std::make_unique<TH1F>("lead_b_pt_ratio_hist", "leading b pair pt ratio", nbins, 0, 6);
    auto lead_l_pt_ratio_hist = std::make_unique<TH1F>("lead_l_pt_ratio_hist", "leading light pair pt ratio", nbins, 0, 6);
    auto sublead_b_pt_ratio_hist = std::make_unique<TH1F>("sublead_b_pt_ratio_hist", "subleading b pair pt ratio", nbins, 0, 6);
    auto sublead_l_pt_ratio_hist = std::make_unique<TH1F>("sublead_l_pt_ratio_hist", "subleading light pair pt ratio", nbins, 0, 6);

    auto w_from_quarks = std::make_unique<TH1F>("w_from_quarks", "W mass from quarks", nbins, 0, 120);
    auto w_from_jets = std::make_unique<TH1F>("w_from_jets", "W mass from jets", nbins/2, 0, 150);

    auto h_from_quarks = std::make_unique<TH1F>("h_from_quarks", "Higgs mass from quarks", nbins, 120, 130);
    auto h_from_jets = std::make_unique<TH1F>("h_from_jets", "Higgs mass from jets", nbins/2, 50, 175);

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

                std::vector<int> partons{ b_idx, bbar_idx, q1_idx, q2_idx };
                
                int b_match = Match(b_idx, genpart, genjet);
                int bbar_match = Match(bbar_idx, genpart, genjet);
                int q1_match = Match(q1_idx, genpart, genjet);
                int q2_match = Match(q2_idx, genpart, genjet);

                // skip event if matching failed
                std::vector<int> matches{b_match, bbar_match, q1_match, q2_match};
                if(std::any_of(matches.begin(), matches.end(), [](int x) { return x == -1;}))
                {
                    ++matching_fails;
                    continue;
                } 
                else
                {
                    ++fill_cnt;
                    double leading_b_pt_ratio = (GenPart_pt[b_idx] > GenPart_pt[bbar_idx]) ? (GenJet_pt[b_match]/GenPart_pt[b_idx]) : (GenJet_pt[bbar_match]/GenPart_pt[bbar_idx]);
                    double leading_l_pt_ratio = (GenPart_pt[q1_idx] > GenPart_pt[q2_idx]) ? (GenJet_pt[q1_match]/GenPart_pt[q1_idx]) : (GenJet_pt[q2_match]/GenPart_pt[q2_idx]);

                    double subleading_b_pt_ratio = (GenPart_pt[b_idx] < GenPart_pt[bbar_idx]) ? (GenJet_pt[b_match]/GenPart_pt[b_idx]) : (GenJet_pt[bbar_match]/GenPart_pt[bbar_idx]);
                    double subleading_l_pt_ratio = (GenPart_pt[q1_idx] < GenPart_pt[q2_idx]) ? (GenJet_pt[q1_match]/GenPart_pt[q1_idx]) : (GenJet_pt[q2_match]/GenPart_pt[q2_idx]);

                    lead_b_pt_ratio_hist->Fill(leading_b_pt_ratio);
                    lead_l_pt_ratio_hist->Fill(leading_l_pt_ratio);

                    sublead_b_pt_ratio_hist->Fill(subleading_b_pt_ratio);
                    sublead_l_pt_ratio_hist->Fill(subleading_l_pt_ratio);

                    std::vector<double> ratios{ leading_b_pt_ratio, subleading_b_pt_ratio, leading_l_pt_ratio, subleading_l_pt_ratio };
                    if (std::any_of(ratios.begin(), ratios.end(), [](double x) { return x > MATCH_PT_THRESH; }))
                    {
                        auto h = EnergyMap(i, genpart, GenPart_genPartIdxMother);
                        auto g1 = Parton(genpart, b_idx);
                        auto g2 = Parton(genpart, bbar_idx);
                        auto g3 = Parton(genpart, q1_idx);
                        auto g4 = Parton(genpart, q2_idx);

                        auto&& graphs = ConeGraphs(genjet, matches, 3);

                        auto c1 = std::make_unique<TCanvas>("c1", "c1");
                        c1->SetGrid();
                        c1->SetTickx();
                        c1->SetTicky();
                        c1->SetLeftMargin(0.15);
                        c1->SetRightMargin(0.15);

                        h->SetStats(0);
                        h->GetXaxis()->SetTitle("phi");
                        h->GetYaxis()->SetTitle("eta");
                        h->Draw("colz");

                        g1->SetLineWidth(2);
                        g1->SetLineColor(kMagenta);
                        g2->SetLineWidth(2);
                        g2->SetLineColor(kMagenta);
                        g3->SetLineWidth(2);
                        g3->SetLineColor(kRed);
                        g4->SetLineWidth(2);
                        g4->SetLineColor(kRed);

                        g1->Draw("same");
                        g2->Draw("same");
                        g3->Draw("same");
                        g4->Draw("same");

                        graphs[0]->Draw("same");
                        graphs[1]->Draw("same");
                        graphs[2]->Draw("same");
                        graphs[3]->Draw("same");
                        
                        // auto leg = std::make_unique<TLegend>(0.15, 0.1, 0.35, 0.3);
                        // leg->AddEntry(g1.get(), "b quarks");
                        // leg->AddEntry(g3.get(), "light quarks");
                        // leg->Draw();

                        c1->SaveAs(Form("Event_%d.png", i));
                        // save_2d_dist(h.get(), Form("Evt_%d_EnergyMap", i), "phi", "eta");
                    }

                    TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                    TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);
                    TLorentzVector bj_p4 = GetP4(genjet, b_match);
                    TLorentzVector bbarj_p4 = GetP4(genjet, bbar_match);

                    TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                    TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                    TLorentzVector lj1_p4 = GetP4(genjet, q1_match);
                    TLorentzVector lj2_p4 = GetP4(genjet, q2_match);

                    w_from_quarks->Fill((lq1_p4 + lq2_p4).M());
                    w_from_jets->Fill((lj1_p4 + lj2_p4).M());
                    h_from_quarks->Fill((bq_p4 + bbarq_p4).M());
                    h_from_jets->Fill((bj_p4 + bbarj_p4).M());
                }
            }
        }
    }

    save_1d_dist(lead_b_pt_ratio_hist.get(), "leading_b_pt_ratio", "jet pt/quark pt");   
    save_1d_dist(lead_l_pt_ratio_hist.get(), "leading_l_pt_ratio", "jet pt/quark pt");

    save_1d_dist(sublead_b_pt_ratio_hist.get(), "subleading_b_pt_ratio", "jet pt/quark pt");   
    save_1d_dist(sublead_l_pt_ratio_hist.get(), "subleading_l_pt_ratio", "jet pt/quark pt");  

    save_1d_dist(w_from_quarks.get(), "w_from_quarks", "[GeV]");   
    save_1d_dist(w_from_jets.get(), "w_from_jets", "[GeV]");  

    save_1d_dist(h_from_quarks.get(), "h_from_quarks", "[GeV]");   
    save_1d_dist(h_from_jets.get(), "h_from_jets", "[GeV]");

    save_1d_stack({w_from_quarks.get(), w_from_jets.get()}, {"quarks", "jets"}, "w_mass", "Hadronically decaying W boson mass", "[GeV]");
    save_1d_stack({h_from_quarks.get(), h_from_jets.get()}, {"quarks", "jets"}, "h_mass", "Higgs boson mass", "[GeV]");     

    std::cout << "nEvents =  " << nEvents << "\n" 
              << "\tnon_empty_sig = " << non_empty_sig << "\n" 
              << "\tvalid_sig = " << valid_sig << "\n"
              << "\tmatching_fails = " << matching_fails << "\n"
              << "\tfill_cnt = " << fill_cnt << "\n";
    // #ifdef DEBUG
    // #endif

    return 0;
}