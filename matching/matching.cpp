#include <iostream>
#include <algorithm>
#include <memory>
#include <iomanip>

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

static constexpr double MATCH_PT_THRESH = 1.2;
static constexpr double DIJET_MASS_THRESH = 95.0;

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

    int overlap_jets = 0;
    int j1lep_overlap = 0;
    int j2lep_overlap = 0;
    int j1nu_overlap = 0;
    int j2nu_overlap = 0;
    int jjlep_overlap = 0;
    int jjnu_overlap = 0;

    int lj1_over_pt_thresh = 0;
    int lj2_over_pt_thresh = 0;
    int dijet_over_mass_thresh = 0;
    int bj_over_pt_thresh = 0;
    int bbarj_over_pt_thresh = 0;

    int failed_parton_cut = 0;

    // hists
    HistManager hm;

    int nbins = 101;
    std::pair<double, double> pt_ratio_range{0.0, 6.0};
    std::pair<double, double> dr_range{0.0, 6.0};
    std::pair<std::string, std::string> mass_labels{"[GeV]", "AU"};
    std::pair<std::string, std::string> pt_ratio_label_labels{"quark pt/jet pt", "AU"};
    std::pair<std::string, std::string> dr_labels{"dR", "AU"};

    // 1D histograms
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
    hm.Add("heavy_higgs_partons", "Heavy higgs mass from partons", mass_labels, {0, 800}, nbins);
    hm.Add("heavy_higgs_jets", "Heavy higgs mass from jets with lepton subtraction", mass_labels, {0, 800}, nbins);
    hm.Add("heavy_higgs_jets_no_corr", "Heavy higgs mass from jets without lepton subtraction", mass_labels, {0, 800}, nbins);
    hm.Add("jet2_jet1_dR", "dR(j2, j1)", dr_labels, dr_range, 15);
    hm.Add("jet2_lep_dR", "dR(j2, lep)", dr_labels, dr_range, 15);
    hm.Add("jet2_nu_dR", "dR(j2, nu)", dr_labels, dr_range, 15);
    hm.Add("jet2_bj_dR", "dR(j2, bj)", dr_labels, dr_range, 15);
    hm.Add("jet2_bbarj_dR", "dR(j2, bbarj)", dr_labels, dr_range, 15);
    hm.Add("jet1_lep_dR", "dR(j1, lep)", dr_labels, dr_range, 15);
    hm.Add("jet1_nu_dR", "dR(j1, nu)", dr_labels, dr_range, 15);
    hm.Add("jet1_bj_dR", "dR(j1, bj)", dr_labels, dr_range, 15);
    hm.Add("jet1_bbarj_dR", "dR(j1, bbarj)", dr_labels, dr_range, 15);
    hm.Add("lead_jet_lep_dR", "dR(lead_jet, lep)", dr_labels, dr_range, 15);
    hm.Add("lead_jet_nu_dR", "dR(lead_jet, nu)", dr_labels, dr_range, 15);
    hm.Add("sublead_jet_lep_dR", "dR(sublead_jet, lep)", dr_labels, dr_range, 15);
    hm.Add("sublead_jet_nu_dR", "dR(sublead_jet, nu)", dr_labels, dr_range, 15);

    std::string sublead_jet_lep_dR("sublead_jet_lep_dR");
    std::string sublead_jet_nu_dR("sublead_jet_nu_dR");
    std::string lead_jet_lep_dR("lead_jet_lep_dR");
    std::string lead_jet_nu_dR("lead_jet_nu_dR");
    std::string jet2_jet1_dR("jet2_jet1_dR");
    std::string jet2_lep_dR("jet2_lep_dR");
    std::string jet2_nu_dR("jet2_nu_dR");
    std::string jet2_bj_dR("jet2_bj_dR");
    std::string jet2_bbarj_dR("jet2_bbarj_dR");
    std::string jet1_lep_dR("jet1_lep_dR");
    std::string jet1_nu_dR("jet1_nu_dR");
    std::string jet1_bj_dR("jet1_bj_dR");
    std::string jet1_bbarj_dR("jet1_bbarj_dR");
    std::string heavy_higgs_partons("heavy_higgs_partons");
    std::string heavy_higgs_jets("heavy_higgs_jets");
    std::string heavy_higgs_jets_no_corr("heavy_higgs_jets_no_corr");
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

    // 2D histograms
    std::pair<double, double> pt_range{0, 250};
    std::pair<int, int> bins{nbins, nbins};
    std::pair<std::string, std::string> pt_label{"quark pt, [GeV]", "jet pt, [GeV]"};
    std::pair<std::string, std::string> dR_pt_label{"dR", "quark_pt/jet_pt"};

    std::string pt_cmp_1("pt_cmp_1");
    std::string pt_cmp_2("pt_cmp_2");
    std::string pt_cmp_lead("pt_cmp_lead");
    std::string pt_cmp_sub("pt_cmp_sub");
    std::string dR_vs_pt_lead_lep("dR_vs_pt_lead_lep");
    std::string dR_vs_pt_lead_nu("dR_vs_pt_lead_nu");
    std::string dR_vs_pt_sublead_lep("dR_vs_pt_sublead_lep");
    std::string dR_vs_pt_sublead_nu("dR_vs_pt_sublead_nu");
    std::string dR_vs_pt_sublead_lead("dR_vs_pt_sublead_lead");
    std::string dR_vs_pt_lead_lead("dR_vs_pt_lead_lead");

    hm.Add(pt_cmp_1, "quark 1 pt vs jet 1 pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_2, "quark 2 pt vs jet 2 pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_lead, "Leading quark pt vs leading jet pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_sub, "Subleading quark pt vs subleading jet pt", pt_label, pt_range, pt_range, bins);
    hm.Add(dR_vs_pt_lead_lep, "dR(lead_jet, lep) vs lead_quark_pt/lead_jet_pt", dR_pt_label, dr_range, dr_range, bins);
    hm.Add(dR_vs_pt_lead_nu, "dR(lead_jet, nu) vs lead_quark_pt/lead_jet_p", dR_pt_label, dr_range, dr_range, bins);
    hm.Add(dR_vs_pt_sublead_lep, "dR(sublead_jet, lep) vs sublead_quark_pt/sublead_jet_pt", dR_pt_label, dr_range, dr_range, bins);
    hm.Add(dR_vs_pt_sublead_nu, "dR(sublead_jet, nu) vs sublead_quark_pt/sublead_jet_p", dR_pt_label, dr_range, dr_range, bins);
    hm.Add(dR_vs_pt_sublead_lead, "dR(sublead_jet, lead_jet) vs sublead_quark_pt/sublead_jet_p", dR_pt_label, dr_range, dr_range, bins);
    hm.Add(dR_vs_pt_lead_lead, "dR(lead_jet, sublead_jet) vs lead_quark_pt/lead_jet_p", dR_pt_label, dr_range, dr_range, bins);

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
                // need to sort matches before doing unique
                // unique removes all except the first CONSECUTIVE duplicates 
                // (see https://en.cppreference.com/w/cpp/algorithm/unique)
                std::sort(matches.begin(), matches.end());
                bool negative_match = std::any_of(matches.begin(), matches.end(), [](int x) { return x == -1;});
                bool same_match = (std::unique(matches.begin(), matches.end()) != matches.end());
                if(negative_match || same_match)
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

                    TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                    TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                    hm.Fill(jet1_nu_dR, lj1_p4.DeltaR(nu_p4));
                    hm.Fill(jet1_lep_dR, lj1_p4.DeltaR(l_p4));
                    hm.Fill(jet1_bj_dR, lj1_p4.DeltaR(bj_p4));
                    hm.Fill(jet1_bbarj_dR, lj1_p4.DeltaR(bbarj_p4));

                    hm.Fill(jet2_jet1_dR, lj2_p4.DeltaR(lj1_p4));
                    hm.Fill(jet2_nu_dR, lj2_p4.DeltaR(nu_p4));
                    hm.Fill(jet2_lep_dR, lj2_p4.DeltaR(l_p4));
                    hm.Fill(jet2_bj_dR, lj2_p4.DeltaR(bj_p4));
                    hm.Fill(jet2_bbarj_dR, lj2_p4.DeltaR(bbarj_p4));

                    TLorentzVector leading_lj = (lj1_p4.Pt() > lj2_p4.Pt()) ? lj1_p4 : lj2_p4;
                    TLorentzVector subleading_lj = (lj1_p4.Pt() < lj2_p4.Pt()) ? lj1_p4 : lj2_p4; 
                    TLorentzVector leading_lq = (lq1_p4.Pt() > lq2_p4.Pt()) ? lq1_p4 : lq2_p4;
                    TLorentzVector subleading_lq = (lq1_p4.Pt() < lq2_p4.Pt()) ? lq1_p4 : lq2_p4; 

                    hm.Fill(pt_cmp_lead, leading_lq.Pt(), leading_lj.Pt());
                    hm.Fill(pt_cmp_sub, subleading_lq.Pt(), subleading_lj.Pt());

                    double lead_lep_dR = leading_lj.DeltaR(l_p4);
                    double lead_nu_dR = leading_lj.DeltaR(nu_p4);
                    double sub_lep_dR = subleading_lj.DeltaR(l_p4);
                    double sub_nu_dR = subleading_lj.DeltaR(nu_p4);

                    hm.Fill(lead_jet_lep_dR, lead_lep_dR);
                    hm.Fill(lead_jet_nu_dR, lead_nu_dR);
                    hm.Fill(sublead_jet_lep_dR, sub_lep_dR);
                    hm.Fill(sublead_jet_nu_dR, sub_nu_dR);

                    hm.Fill(dR_vs_pt_lead_lep, lead_lep_dR, leading_lq.Pt()/leading_lj.Pt());
                    hm.Fill(dR_vs_pt_lead_nu, lead_nu_dR, leading_lq.Pt()/leading_lj.Pt());
                    hm.Fill(dR_vs_pt_sublead_lep, sub_lep_dR, subleading_lq.Pt()/subleading_lj.Pt());
                    hm.Fill(dR_vs_pt_sublead_nu, sub_nu_dR, subleading_lq.Pt()/subleading_lj.Pt());

                    hm.Fill(dR_vs_pt_sublead_lead, subleading_lj.DeltaR(leading_lj), subleading_lq.Pt()/subleading_lj.Pt());
                    hm.Fill(dR_vs_pt_lead_lead, leading_lj.DeltaR(subleading_lj), leading_lq.Pt()/leading_lj.Pt());

                    // if (lq1_p4.Pt() < 20 || lq2_p4.Pt() < 20)
                    // {
                    //     ++failed_parton_cut;
                    //     continue;
                    // }

                    double hh_mass_jets_no_corr = (lj1_p4 + lj2_p4 + l_p4 + nu_p4 + bj_p4 + bbarj_p4).M();
                    hm.Fill(heavy_higgs_jets_no_corr, hh_mass_jets_no_corr);

                    bool overlapping_jets = lj1_p4.DeltaR(lj2_p4) < 2*DR_THRESH;
                    if (overlapping_jets)
                    {
                        ++overlap_jets;
                        hm.Fill(overlap_pt_ratio_1, lq1_p4.Pt()/lj1_p4.Pt());
                        hm.Fill(overlap_pt_ratio_2, lq2_p4.Pt()/lj2_p4.Pt());
                    }

                    // subtract lepton from jet if they overlap
                    double dR_lj1_lep = lj1_p4.DeltaR(l_p4);
                    double dR_lj2_lep = lj2_p4.DeltaR(l_p4);
                    if (dR_lj1_lep < DR_THRESH)
                    {
                        ++j1lep_overlap;
                        if (lj1_p4.E() > l_p4.E()) lj1_p4 -= l_p4;
                    }
                    if (dR_lj2_lep < DR_THRESH) 
                    {
                        ++j2lep_overlap;
                        if (lj2_p4.E() > l_p4.E()) lj2_p4 -= l_p4;
                    }

                    // subtract neutrino from jet if they overlap
                    double dR_lj1_nu = lj1_p4.DeltaR(nu_p4);
                    double dR_lj2_nu = lj2_p4.DeltaR(nu_p4);
                    if (dR_lj1_nu < DR_THRESH)
                    {
                        ++j1nu_overlap;
                        // if (lj1_p4.E() > nu_p4.E()) lj1_p4 -= nu_p4;
                    }
                    if (dR_lj2_nu < DR_THRESH) 
                    {
                        ++j2nu_overlap;
                        // if (lj2_p4.E() > nu_p4.E()) lj2_p4 -= nu_p4;;
                    }

                    if (dR_lj1_lep < DR_THRESH && dR_lj2_lep < DR_THRESH) ++jjlep_overlap;
                    if (dR_lj1_nu < DR_THRESH && dR_lj2_nu < DR_THRESH) ++jjnu_overlap;

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

                    lj1_over_pt_thresh += (lj1_p4.Pt()/lq1_p4.Pt() > MATCH_PT_THRESH) ? 1 : 0;
                    lj2_over_pt_thresh += (lj2_p4.Pt()/lq2_p4.Pt() > MATCH_PT_THRESH) ? 1 : 0;
                    bj_over_pt_thresh += (bj_p4.Pt()/bq_p4.Pt() > MATCH_PT_THRESH) ? 1 : 0;
                    bbarj_over_pt_thresh += (bbarj_p4.Pt()/bbarq_p4.Pt() > MATCH_PT_THRESH) ? 1 : 0;

                    hm.Fill(pt_cmp_1, lq1_p4.Pt(), lj1_p4.Pt());
                    hm.Fill(pt_cmp_2, lq2_p4.Pt(), lj2_p4.Pt());

                    hm.Fill(w_from_quarks, (lq1_p4 + lq2_p4).M());
                    hm.Fill(w_from_jets, (lj1_p4 + lj2_p4).M());
                    hm.Fill(h_from_quarks, (bq_p4 + bbarq_p4).M());
                    hm.Fill(h_from_jets, (bj_p4 + bbarj_p4).M());

                    double dijet_mass = (lj1_p4 + lj2_p4).M();
                    dijet_over_mass_thresh += (dijet_mass > DIJET_MASS_THRESH) ? 1 : 0;

                    double hh_mass_part = (lq1_p4 + lq2_p4 + l_p4 + nu_p4 + bq_p4 + bbarq_p4).M();
                    hm.Fill(heavy_higgs_partons, hh_mass_part);
                    double hh_mass_jets = (lj1_p4 + lj2_p4 + l_p4 + nu_p4 + bj_p4 + bbarj_p4).M();
                    hm.Fill(heavy_higgs_jets, hh_mass_jets);
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
    hm.DrawStack({heavy_higgs_partons, heavy_higgs_jets, heavy_higgs_jets_no_corr}, "Heavy higgs mass", "hh_mass");
    hm.DrawStack({heavy_higgs_jets, heavy_higgs_jets_no_corr}, "Lepton subtraction effect on hh mass", "hh_mass_jets");
    hm.DrawStack({jet2_jet1_dR, jet2_lep_dR, jet2_nu_dR, jet2_bj_dR, jet2_bbarj_dR}, "dR(lj2, *)", "jet2_dR");
    hm.DrawStack({jet2_jet1_dR, jet1_lep_dR, jet1_nu_dR, jet1_bj_dR, jet1_bbarj_dR}, "dR(lj1, *)", "jet1_dR");
    hm.DrawStack({lead_jet_lep_dR, lead_jet_nu_dR, sublead_jet_lep_dR, sublead_jet_nu_dR}, "dR between jets and leptons", "jets_dR");

    std::cout << std::setprecision(3);
    std::cout << "nEvents = " << nEvents << "\n" 
              << "\tnon_empty_sig = " << non_empty_sig << "\n" 
              << "\tvalid_sig = " << valid_sig << "\n"
              << "\tmatching_fails = " << matching_fails << "\n"
              << "\taccep_evts = " << accep_evts << "\n"
              << "\tfailed_parton_cut = " << failed_parton_cut << "\n"
              << "\tlj1_over_pt_thresh = " << lj1_over_pt_thresh << " (" << 1.0*lj1_over_pt_thresh/accep_evts*100 << "%)" << "\n"
              << "\tlj2_over_pt_thresh = " << lj2_over_pt_thresh << " (" << 1.0*lj2_over_pt_thresh/accep_evts*100 << "%)" << "\n"
              << "\tbj_over_pt_thresh = " << bj_over_pt_thresh << " (" << 1.0*bj_over_pt_thresh/accep_evts*100 << "%)" << "\n"
              << "\tbbarj_over_pt_thresh = " << bbarj_over_pt_thresh << " (" << 1.0*bbarj_over_pt_thresh/accep_evts*100 << "%)" << "\n"
              << "\toverlap_jets = " << overlap_jets << " (" << 1.0*overlap_jets/accep_evts*100 << "%)" << "\n"
              << "\tj1lep_overlap = " << j1lep_overlap << " (" << 1.0*j1lep_overlap/accep_evts*100 << "%)" << "\n"
              << "\tj2lep_overlap = " << j2lep_overlap << " (" << 1.0*j2lep_overlap/accep_evts*100 << "%)" << "\n"
              << "\tj1nu_overlap = " << j1nu_overlap << " (" << 1.0*j1nu_overlap/accep_evts*100 << "%)" << "\n"
              << "\tj2nu_overlap = " << j2nu_overlap << " (" << 1.0*j2nu_overlap/accep_evts*100 << "%)" << "\n"
              << "\tjjlep_overlap = " << jjlep_overlap << "\n"
              << "\tjjnu_overlap = " << jjnu_overlap << "\n"
              << "\tdijet_over_mass_thresh = " << dijet_over_mass_thresh << " (" << 1.0*dijet_over_mass_thresh/accep_evts*100 << "%)" << "\n";

    return 0;
}