#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#ifdef DEBUG
    #include "TError.h"
#endif

#include "MatchingTools.hpp"
#include "HistManager.hpp"

static constexpr int MAX_GENJET = 21;
static constexpr int MAX_GENPART = 270;
static constexpr double UP_PT_THRESH = 1.2;
static constexpr double LOW_PT_THRESH = 0.8;
static constexpr double DIJET_MASS_THRESH = 95.0;
static constexpr double MIN_PARTON_PT = 25.0;

int main()
{
    #ifdef DEBUG
        gErrorIgnoreLevel = 3000;
    #endif
    
    std::cout << std::setprecision(3);
    auto start = std::chrono::system_clock::now();
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    // TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-450.root");
    // TFile *myFile = TFile::Open("NanoAOD_600GeV_1000Events.root");
    // TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-600_narrow_13TeV-madgraph.root");
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
    int matched_events = 0;
    int accepted_events = 0;

    int overlap_jets = 0;
    int j1lep_overlap = 0;
    int j2lep_overlap = 0;

    int lj1_above_pt_thresh = 0;
    int lj2_above_pt_thresh = 0;
    int dijet_above_mass_thresh = 0;
    int bj_above_pt_thresh = 0;
    int bbarj_above_pt_thresh = 0;
    int lj1_under_pt_thresh = 0;
    int lj2_under_pt_thresh = 0;
    int bj_under_pt_thresh = 0;
    int bbarj_under_pt_thresh = 0;

    int failed_min_dR = 0;
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
    std::string response_all("response_all");
    std::string resolution_all("resolution_all");
    std::string heavy_higgs_partons("heavy_higgs_partons");
    std::string heavy_higgs_jets("heavy_higgs_jets");
    std::string heavy_higgs_jets_wnu("heavy_higgs_jets_wnu");
    std::string w_from_quarks("w_from_quarks");
    std::string w_from_jets("w_from_jets");
    std::string h_from_quarks("h_from_quarks");
    std::string h_from_jets("h_from_jets");
    std::string h_from_jets_wnu("h_from_jets_wnu");

    hm.Add(w_from_quarks, "W mass from quarks", mass_labels, {0, 120}, nbins);
    hm.Add(w_from_jets, "W mass from jets", mass_labels, {0, 150}, nbins);
    hm.Add(h_from_quarks, "Higgs mass from quarks", mass_labels, {120, 130}, nbins);
    hm.Add(h_from_jets, "Higgs mass from jets", mass_labels, {0, 200}, nbins);
    hm.Add(h_from_jets_wnu, "Higgs mass from jets with neutrinos", mass_labels, {0, 200}, nbins);
    hm.Add(heavy_higgs_partons, "Heavy higgs mass from partons", mass_labels, {0, 800}, nbins);
    hm.Add(heavy_higgs_jets, "Heavy higgs mass from jets with lepton subtraction", mass_labels, {0, 800}, nbins);
    hm.Add(heavy_higgs_jets_wnu, "Heavy higgs mass from jets with lepton subtraction and missing neutrinos added", mass_labels, {0, 800}, nbins);
    hm.Add(response_all, "Response", {"jet pt/quark pt", "AU"}, {0, 2}, nbins);
    hm.Add(resolution_all, "Resolution", {"(jet pt - quark pt)/quark pt", "AU"}, {-1, 1}, nbins);


    // 2D histograms
    std::pair<double, double> pt_range{0, 350};
    std::pair<int, int> bins{nbins, nbins};
    std::pair<std::string, std::string> pt_label{"quark pt, [GeV]", "jet pt, [GeV]"};
    std::pair<std::string, std::string> dR_ptr_label{"dR", "quark_pt/jet_pt"};
    std::pair<std::string, std::string> dR_pt_label{"dR", "quark_pt"};

    std::string pt_cmp_l1("pt_cmp_l1");
    std::string pt_cmp_l2("pt_cmp_l2");
    std::string pt_cmp_b("pt_cmp_b");
    std::string pt_cmp_bbar("pt_cmp_bbar");

    hm.Add(pt_cmp_l1, "quark 1 pt vs jet 1 pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_l2, "quark 2 pt vs jet 2 pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_b, "b quark pt vs b jet pt", pt_label, pt_range, pt_range, bins);
    hm.Add(pt_cmp_bbar, "bbar quark pt vs bar jet pt", pt_label, pt_range, pt_range, bins);

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
                    ++matched_events;

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

                    std::vector<TLorentzVector> partons{ lq1_p4, lq2_p4, bq_p4, bbarq_p4 };
                    if (std::any_of(partons.begin(), partons.end(), [](TLorentzVector const& p) { return p.Pt() < MIN_PARTON_PT; }))
                    {
                        ++failed_parton_cut;
                        // continue;
                    }

                    std::vector<TLorentzVector> parts{lj1_p4, lj2_p4, bj_p4, bbarj_p4, l_p4, nu_p4};
                    if (MinDeltaR(parts) < DR_THRESH)
                    {
                        ++failed_min_dR;
                        // continue;
                    }

                    ++accepted_events;

                    #ifdef DEBUG
                        TLorentzVector& leading_lj = (lj1_p4.Pt() > lj2_p4.Pt()) ? lj1_p4 : lj2_p4;
                        TLorentzVector& subleading_lj = (lj1_p4.Pt() < lj2_p4.Pt()) ? lj1_p4 : lj2_p4; 
                        TLorentzVector& leading_lq = (lq1_p4.Pt() > lq2_p4.Pt()) ? lq1_p4 : lq2_p4;
                        TLorentzVector& subleading_lq = (lq1_p4.Pt() < lq2_p4.Pt()) ? lq1_p4 : lq2_p4; 

                        TLorentzVector& leading_bj = (bj_p4.Pt() > bbarj_p4.Pt()) ? bj_p4 : bbarj_p4;
                        TLorentzVector& subleading_bj = (bj_p4.Pt() < bbarj_p4.Pt()) ? bj_p4 : bbarj_p4;
                        TLorentzVector& leading_bq = (bq_p4.Pt() > bbarq_p4.Pt()) ? bq_p4 : bbarq_p4;
                        TLorentzVector& subleading_bq = (bq_p4.Pt() < bbarq_p4.Pt()) ? bq_p4 : bbarq_p4;

                        double leading_b_pt_ratio = leading_bj.Pt()/leading_bq.Pt();
                        double leading_l_pt_ratio = leading_lj.Pt()/leading_lq.Pt();
                        double subleading_b_pt_ratio = subleading_bj.Pt()/subleading_bq.Pt();
                        double subleading_l_pt_ratio = subleading_lj.Pt()/subleading_lq.Pt();

                        std::vector<double> ratios{ leading_b_pt_ratio, subleading_b_pt_ratio, leading_l_pt_ratio, subleading_l_pt_ratio };
                        if (std::any_of(ratios.begin(), ratios.end(), [](double x) { return (x > UP_PT_THRESH || x < LOW_PT_THRESH); }))
                        {
                            MatchIndex mi = {{b_idx, b_match}, {bbar_idx, bbar_match}, {q1_idx, q1_match}, {q2_idx, q2_match}};
                            MatchKinematics mk = {genpart, genjet};
                            std::pair<int*, int*> ptrs = {GenPart_genPartIdxMother, GenPart_pdgId};
                            DrawEventMap(mk, mi, i, ptrs);
                        }
                    #endif  

                    bool overlapping_jets = lj1_p4.DeltaR(lj2_p4) < 2*DR_THRESH;
                    if (overlapping_jets)
                    {
                        ++overlap_jets;
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

                    hm.Fill(response_all, bj_p4.Pt()/bq_p4.Pt());
                    hm.Fill(response_all, bbarj_p4.Pt()/bbarq_p4.Pt());
                    hm.Fill(response_all, lj1_p4.Pt()/lq1_p4.Pt());
                    hm.Fill(response_all, lj2_p4.Pt()/lq2_p4.Pt());

                    hm.Fill(resolution_all, (bj_p4.Pt() - bq_p4.Pt())/bq_p4.Pt());
                    hm.Fill(resolution_all, (bbarj_p4.Pt() - bbarq_p4.Pt())/bbarq_p4.Pt());
                    hm.Fill(resolution_all, (lj1_p4.Pt() - lq1_p4.Pt())/lq1_p4.Pt());
                    hm.Fill(resolution_all, (lj2_p4.Pt() - lq2_p4.Pt())/lq2_p4.Pt());

                    hm.Fill(w_from_quarks, (lq1_p4 + lq2_p4).M());
                    hm.Fill(w_from_jets, (lj1_p4 + lj2_p4).M());
                    hm.Fill(h_from_quarks, (bq_p4 + bbarq_p4).M());
                    hm.Fill(h_from_jets, (bj_p4 + bbarj_p4).M());

                    double dijet_mass = (lj1_p4 + lj2_p4).M();
                    dijet_above_mass_thresh += (dijet_mass > DIJET_MASS_THRESH) ? 1 : 0;

                    double hh_mass_part = (lq1_p4 + lq2_p4 + l_p4 + nu_p4 + bq_p4 + bbarq_p4).M();
                    hm.Fill(heavy_higgs_partons, hh_mass_part);
                    double hh_mass_jets = (lj1_p4 + lj2_p4 + l_p4 + nu_p4 + bj_p4 + bbarj_p4).M();
                    hm.Fill(heavy_higgs_jets, hh_mass_jets);

                    // determines if particle at index idx is a neutrino
                    auto NotNu = [&GenPart_pdgId](int idx) { return !IsNeutrino(GenPart_pdgId[idx]); };

                    // find neutrinos from b jets
                    auto b_neutrinos = GetStableDescendants(b_idx, GenPart_genPartIdxMother, nGenPart);
                    b_neutrinos.erase(std::remove_if(b_neutrinos.begin(), b_neutrinos.end(), NotNu), b_neutrinos.end());
                    auto bbar_neutrinos = GetStableDescendants(bbar_idx, GenPart_genPartIdxMother, nGenPart);
                    bbar_neutrinos.erase(std::remove_if(bbar_neutrinos.begin(), bbar_neutrinos.end(), NotNu), bbar_neutrinos.end());

                    // 4-momentum generator
                    auto GenP4 = [&genpart](int idx) { return GetP4(genpart, idx); };

                    std::vector<TLorentzVector> nu_b_p4;
                    std::transform(b_neutrinos.begin(), b_neutrinos.end(), std::back_inserter(nu_b_p4), GenP4);
                    TLorentzVector tot_nu_b_p4 = std::accumulate(nu_b_p4.begin(), nu_b_p4.end(), TLorentzVector());
                    bj_p4 += tot_nu_b_p4;

                    std::vector<TLorentzVector> nu_bbar_p4;
                    std::transform(bbar_neutrinos.begin(), bbar_neutrinos.end(), std::back_inserter(nu_bbar_p4), GenP4);
                    TLorentzVector tot_nu_bbar_p4 = std::accumulate(nu_bbar_p4.begin(), nu_bbar_p4.end(), TLorentzVector());
                    bbarj_p4 += tot_nu_bbar_p4;
                    hh_mass_jets = (lj1_p4 + lj2_p4 + l_p4 + nu_p4 + bj_p4 + bbarj_p4).M();
                    hm.Fill(heavy_higgs_jets_wnu, hh_mass_jets);
                    hm.Fill(h_from_jets_wnu, (bj_p4 + bbarj_p4).M());

                    #ifdef DEBUG
                        std::vector<double> pt_ratios{ lj1_p4.Pt()/lq1_p4.Pt(), lj2_p4.Pt()/lq2_p4.Pt(), bj_p4.Pt()/bq_p4.Pt(), bbarj_p4.Pt()/bbarq_p4.Pt() };
                        if (std::any_of(pt_ratios.begin(), pt_ratios.end(), [](double x) { return x > UP_PT_THRESH; }))
                        {
                            std::cout << "Event " << i << ":\n";
                            std::cout << "DeltaR:\n";
                            std::cout << "dR(q1, q2) = " << lq1_p4.DeltaR(lq2_p4) << "\n"
                                    << "dR(q1, bq) = " << lq1_p4.DeltaR(bq_p4) << "\n"
                                    << "dR(q2, bq) = " << lq2_p4.DeltaR(bq_p4) << "\n"
                                    << "dR(q1, bbarq) = " << lq1_p4.DeltaR(bbarq_p4) << "\n"
                                    << "dR(q2, bbarq) = " << lq2_p4.DeltaR(bbarq_p4) << "\n"
                                    << "dR(bq, bbarq) = " << bq_p4.DeltaR(bbarq_p4) << "\n";
                            std::cout << "\n";

                            std::cout << "j1 = "; Print(lj1_p4);
                            std::cout << "q1 = "; Print(lq1_p4);
                            std::cout << "ratio = " << lj1_p4.Pt()/lq1_p4.Pt() << "\n\n";

                            std::cout << "j2 = "; Print(lj2_p4);
                            std::cout << "q2 = "; Print(lq2_p4);
                            std::cout << "ratio = " << lj2_p4.Pt()/lq2_p4.Pt() << "\n\n";

                            std::cout << "bj = "; Print(bj_p4);
                            std::cout << "bq = "; Print(bq_p4);
                            std::cout << "ratio = " << bj_p4.Pt()/bq_p4.Pt() << "\n\n";

                            std::cout << "bbarj = "; Print(bbarj_p4);
                            std::cout << "bbarq = "; Print(bbarq_p4);
                            std::cout << "ratio = " << bbarj_p4.Pt()/bbarq_p4.Pt() << "\n\n";

                            // find stable daughters of all quarks and summ all their momenta
                            auto PrintP4 = [&genpart](int i)
                            {
                                TLorentzVector p4 = GetP4(genpart, i);
                                std::cout << "\t";
                                Print(p4);
                            };

                            TLorentzVector sum;
                            std::vector<int> b_stable_desc = GetStableDescendants(b_idx, GenPart_genPartIdxMother, nGenPart);
                            std::cout << "b: ";
                            std::for_each(b_stable_desc.begin(), b_stable_desc.end(), [&GenPart_pdgId](int i){ std::cout << GenPart_pdgId[i] << " "; });
                            sum = std::transform_reduce(b_stable_desc.begin(), b_stable_desc.end(), TLorentzVector(), std::plus{}, GenP4);
                            std::cout << "\n";
                            std::for_each(b_stable_desc.begin(), b_stable_desc.end(), PrintP4);
                            std::cout << "sum = "; Print(sum);

                            std::vector<int> bbar_stable_desc = GetStableDescendants(bbar_idx, GenPart_genPartIdxMother, nGenPart);
                            std::cout << "bbar: ";
                            std::for_each(bbar_stable_desc.begin(), bbar_stable_desc.end(), [&GenPart_pdgId](int i){ std::cout << GenPart_pdgId[i] << " "; });
                            sum = std::transform_reduce(bbar_stable_desc.begin(), bbar_stable_desc.end(), TLorentzVector(), std::plus{}, GenP4);
                            std::cout << "\n";
                            std::for_each(bbar_stable_desc.begin(), bbar_stable_desc.end(), PrintP4);
                            std::cout << "sum = "; Print(sum);

                            std::vector<int> q1_stable_desc = GetStableDescendants(q1_idx, GenPart_genPartIdxMother, nGenPart);
                            std::cout << "q1: ";
                            std::for_each(q1_stable_desc.begin(), q1_stable_desc.end(), [&GenPart_pdgId](int i){ std::cout << GenPart_pdgId[i] << " "; });
                            sum = std::transform_reduce(q1_stable_desc.begin(), q1_stable_desc.end(), TLorentzVector(), std::plus{}, GenP4);
                            std::cout << "\n";
                            std::for_each(q1_stable_desc.begin(), q1_stable_desc.end(), PrintP4);
                            std::cout << "sum = "; Print(sum);

                            std::vector<int> q2_stable_desc = GetStableDescendants(q2_idx, GenPart_genPartIdxMother, nGenPart);
                            std::cout << "q2: ";
                            std::for_each(q2_stable_desc.begin(), q2_stable_desc.end(), [&GenPart_pdgId](int i){ std::cout << GenPart_pdgId[i] << " "; });
                            sum = std::transform_reduce(q2_stable_desc.begin(), q2_stable_desc.end(), TLorentzVector(), std::plus{}, GenP4);
                            std::cout << "\n";
                            std::for_each(q2_stable_desc.begin(), q2_stable_desc.end(), PrintP4);
                            std::cout << "sum = "; Print(sum);
                            
                            std::cout << "=============================================================================\n";
                        }
                    #endif

                    lj1_above_pt_thresh += (lj1_p4.Pt()/lq1_p4.Pt() > UP_PT_THRESH) ? 1 : 0;
                    lj2_above_pt_thresh += (lj2_p4.Pt()/lq2_p4.Pt() > UP_PT_THRESH) ? 1 : 0;
                    bj_above_pt_thresh += (bj_p4.Pt()/bq_p4.Pt() > UP_PT_THRESH) ? 1 : 0;
                    bbarj_above_pt_thresh += (bbarj_p4.Pt()/bbarq_p4.Pt() > UP_PT_THRESH) ? 1 : 0;

                    lj1_under_pt_thresh += (lj1_p4.Pt()/lq1_p4.Pt() < LOW_PT_THRESH) ? 1 : 0;
                    lj2_under_pt_thresh += (lj2_p4.Pt()/lq2_p4.Pt() < LOW_PT_THRESH) ? 1 : 0;
                    bj_under_pt_thresh += (bj_p4.Pt()/bq_p4.Pt() < LOW_PT_THRESH) ? 1 : 0;
                    bbarj_under_pt_thresh += (bbarj_p4.Pt()/bbarq_p4.Pt() < LOW_PT_THRESH) ? 1 : 0;

                    hm.Fill(pt_cmp_l1, lq1_p4.Pt(), lj1_p4.Pt());
                    hm.Fill(pt_cmp_l2, lq2_p4.Pt(), lj2_p4.Pt());
                    hm.Fill(pt_cmp_b, bq_p4.Pt(), bj_p4.Pt());
                    hm.Fill(pt_cmp_bbar, bbarq_p4.Pt(), bbarj_p4.Pt());
                }
            }
        }
    }

    std::cout << "Finished processing\n";
    hm.Draw();
    hm.DrawStack({w_from_quarks, w_from_jets}, "W mass", "w_mass");
    hm.DrawStack({h_from_jets, h_from_jets_wnu}, "Higgs mass", "higgs_mass");
    hm.DrawStack({heavy_higgs_jets, heavy_higgs_jets_wnu}, "Heavy higgs mass", "hh_mass");

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n" 
              << "\tnon_empty_sig = " << non_empty_sig << "\n" 
              << "\tvalid_sig = " << valid_sig << "\n"
              << "\tmatching_fails = " << matching_fails << "\n"
              << "\tmatched_events = " << matched_events << "\n"
              << "\tfailed_min_dR = " << failed_min_dR << " (" << 1.0*failed_min_dR/matched_events*100 << "%)" << "\n"
              << "\tfailed_parton_cut = " << failed_parton_cut << " (" << 1.0*failed_parton_cut/matched_events*100 << "%)" << "\n"
              << "\taccepted_events = " << accepted_events << " (" << 1.0*accepted_events/matched_events*100 << "%)" << "\n"
              << "\tlj1_above_pt_thresh = " << lj1_above_pt_thresh << " (" << 1.0*lj1_above_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tlj2_above_pt_thresh = " << lj2_above_pt_thresh << " (" << 1.0*lj2_above_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tbj_above_pt_thresh = " << bj_above_pt_thresh << " (" << 1.0*bj_above_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tbbarj_above_pt_thresh = " << bbarj_above_pt_thresh << " (" << 1.0*bbarj_above_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tlj1_under_pt_thresh = " << lj1_under_pt_thresh << " (" << 1.0*lj1_under_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tlj2_under_pt_thresh = " << lj2_under_pt_thresh << " (" << 1.0*lj2_under_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tbj_under_pt_thresh = " << bj_under_pt_thresh << " (" << 1.0*bj_under_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\tbbarj_under_pt_thresh = " << bbarj_under_pt_thresh << " (" << 1.0*bbarj_under_pt_thresh/accepted_events*100 << "%)" << "\n"
              << "\toverlap_jets = " << overlap_jets << " (" << 1.0*overlap_jets/accepted_events*100 << "%)" << "\n"
              << "\tj1lep_overlap = " << j1lep_overlap << " (" << 1.0*j1lep_overlap/accepted_events*100 << "%)" << "\n"
              << "\tj2lep_overlap = " << j2lep_overlap << " (" << 1.0*j2lep_overlap/accepted_events*100 << "%)" << "\n"
              << "\tdijet_above_mass_thresh = " << dijet_above_mass_thresh << " (" << 1.0*dijet_above_mass_thresh/accepted_events*100 << "%)" << "\n";
    return 0;
}