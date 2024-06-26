#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <chrono>
#include <cassert>

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

static constexpr int MAX_AK4_GENJET = 21;
static constexpr int MAX_AK8_GENJET = 7;
static constexpr int MAX_GENPART = 270;

int main()
{
    auto start = std::chrono::system_clock::now();
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    UInt_t          nGenPart;
    Float_t         GenPart_eta[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_mass[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_phi[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_pt[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_pdgId[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_status[MAX_GENPART];   //[nGenPart]

    // ak4 jets
    UInt_t          nGenJetAK4;
    Float_t         GenJetAK4_eta[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_mass[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_phi[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_pt[MAX_AK4_GENJET];   //[nGenJetAK4]
    Int_t           GenJetAK4_partonFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    UChar_t         GenJetAK4_hadronFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    
    // ak8 jets
    // UInt_t          nGenJetAK8;
    // Float_t         GenJetAK8_eta[MAX_AK8_GENJET];   //[nGenJetAK8]
    // Float_t         GenJetAK8_mass[MAX_AK8_GENJET];   //[nGenJetAK8]
    // Float_t         GenJetAK8_phi[MAX_AK8_GENJET];   //[nGenJetAK8]
    // Float_t         GenJetAK8_pt[MAX_AK8_GENJET];   //[nGenJetAK8]
    // Int_t           GenJetAK8_partonFlavour[MAX_AK8_GENJET];   //[nGenJetAK8]
    // UChar_t         GenJetAK8_hadronFlavour[MAX_AK8_GENJET];   //[nGenJetAK8]

    // ak8 subjets
    // UInt_t          nSubGenJetAK8;
    // Float_t         SubGenJetAK8_eta[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    // Float_t         SubGenJetAK8_mass[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    // Float_t         SubGenJetAK8_phi[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    // Float_t         SubGenJetAK8_pt[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]

    Float_t         GenMET_phi;
    Float_t         GenMET_pt;


    // gen particles 
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

    // ak4 gen jets
    myTree->SetBranchAddress("nGenJet", &nGenJetAK4);
    myTree->SetBranchAddress("GenJet_eta", &GenJetAK4_eta);
    myTree->SetBranchAddress("GenJet_mass", &GenJetAK4_mass);
    myTree->SetBranchAddress("GenJet_phi", &GenJetAK4_phi);
    myTree->SetBranchAddress("GenJet_pt", &GenJetAK4_pt);
    myTree->SetBranchAddress("GenJet_partonFlavour", &GenJetAK4_partonFlavour);
    myTree->SetBranchAddress("GenJet_hadronFlavour", &GenJetAK4_hadronFlavour);

    // see https://stackoverflow.com/questions/15585267/cout-not-printing-unsigned-char to use GenJet_hadronFlavour
    // explicitly cast to unsigned
    // also can be fixed by unary +
    // unsigned is understood as unsifned int by compiler

    // ak8 gen jets
    // myTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    // myTree->SetBranchAddress("GenJetAK8_eta", &GenJetAK8_eta);
    // myTree->SetBranchAddress("GenJetAK8_mass", &GenJetAK8_mass);
    // myTree->SetBranchAddress("GenJetAK8_phi", &GenJetAK8_phi);
    // myTree->SetBranchAddress("GenJetAK8_pt", &GenJetAK8_pt);
    // myTree->SetBranchAddress("GenJetAK8_partonFlavour", &GenJetAK8_partonFlavour);
    // myTree->SetBranchAddress("GenJetAK8_hadronFlavour", &GenJetAK8_hadronFlavour);

    // ak8 hen subjets
    // myTree->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    // myTree->SetBranchAddress("SubGenJetAK8_eta", &SubGenJetAK8_eta);
    // myTree->SetBranchAddress("SubGenJetAK8_mass", &SubGenJetAK8_mass);
    // myTree->SetBranchAddress("SubGenJetAK8_phi", &SubGenJetAK8_phi);
    // myTree->SetBranchAddress("SubGenJetAK8_pt", &SubGenJetAK8_pt);

    myTree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    myTree->SetBranchAddress("GenMET_pt", &GenMET_pt);

    int nEvents = myTree->GetEntries();

    // counters

    int has_bbWW_decay = 0;
    int has_tau = 0;
    int has_only_emu = 0;
    int has_at_least_4_primary_jets = 0;
    int has_at_least_4_clean_jets = 0;

    int has_at_least_1_bflav_jet = 0;
    int has_at_least_2_bflav_jet = 0;
    int has_at_least_2_bflav_jet_passing_pt = 0;
    int has_at_least_2_bflav_jet_passing_eta = 0;
    int has_accept_lep = 0;

    int matched_events = 0;
    int inconsistent_b_jet_match = 0;
    int inconsistent_light_jet_match = 0;
    int accepted_events = 0;

    // hists
    HistManager hm;

    // 1D histograms
    std::string higgs_mass_jets("higgs_mass_jets");
    std::string higgs_mass_jets_wnu("higgs_mass_jets_wnu");
    std::string w_mass_jets("w_mass_jets");
    std::string w_mass_jets_wnu("w_mass_jets_wnu");
    std::string pt_ratio_b("pt_ratio_b");
    std::string pt_ratio_bbar("pt_ratio_bbar");
    std::string pt_ratio_b_wnu("pt_ratio_b_wnu");
    std::string pt_ratio_bbar_wnu("pt_ratio_bbar_wnu");
    std::string pt_ratio_1("pt_ratio_1");
    std::string pt_ratio_2("pt_ratio_2");
    std::string pt_ratio_1_wnu("pt_ratio_1_wnu");
    std::string pt_ratio_2_wnu("pt_ratio_2_wnu");

    hm.Add(higgs_mass_jets, "Higgs mass", {"Higgs mass, [GeV]", "Count"}, {40, 160}, 40);
    hm.Add(higgs_mass_jets_wnu, "Higgs mass", {"Higgs mass, [GeV]", "Count"}, {60, 160}, 40);
    hm.Add(w_mass_jets, "W mass", {"W mass, [GeV]", "Count"}, {0, 120}, 40);
    hm.Add(w_mass_jets_wnu, "W mass", {"W mass, [GeV]", "Count"}, {0, 120}, 40);
    hm.Add(pt_ratio_b, "(Jet_pt - quark_pt)/quark_pt: b", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_bbar, "(Jet_pt - quark_pt)/quark_pt: bbar", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_b_wnu, "(Jet_pt - quark_pt)/quark_pt: b (neutrino added)", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_bbar_wnu, "(Jet_pt - quark_pt)/quark_pt: bbar (neutrino added)", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_1, "(Jet_pt - quark_pt)/quark_pt: 1", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_2, "(Jet_pt - quark_pt)/quark_pt: 2", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_1_wnu, "(Jet_pt - quark_pt)/quark_pt: 1 (neutrino added)", {"ratio", "Count"}, {-3, 3}, 50);
    hm.Add(pt_ratio_2_wnu, "(Jet_pt - quark_pt)/quark_pt: 2 (neutrino added)", {"ratio", "Count"}, {-3, 3}, 50);

    // 2D histograms
    std::string met_vs_all_nu("met_vs_all_nu");
    std::string met_vs_all_nu_px("met_vs_all_nu_px");
    std::string met_vs_all_nu_py("met_vs_all_nu_py");
    std::string jet_vs_quark_b("jet_vs_quark_b");
    std::string jet_vs_quark_bbar("jet_vs_quark_bbar");
    std::string jet_vs_quark_b_wnu("jet_vs_quark_b_wnu");
    std::string jet_vs_quark_bbar_wnu("jet_vs_quark_bbar_wnu");
    std::string jet_vs_quark_q1("jet_vs_quark_q1");
    std::string jet_vs_quark_q2("jet_vs_quark_q2");
    std::string jet_vs_quark_q1_wnu("jet_vs_quark_q1_wnu");
    std::string jet_vs_quark_q2_wnu("jet_vs_quark_q2_wnu");
    std::string hadW_vs_lepW("hadW_vs_lepW");

    hm.Add(met_vs_all_nu, "MET vs sum of all neutrinos in the event", {"MET pt, [GeV]", "all nu pt, [GeV]"}, {0, 250}, {0, 250}, {50, 50});
    hm.Add(met_vs_all_nu_px, "MET vs sum of all neutrinos in the event", {"MET px, [GeV]", "all nu px, [GeV]"}, {0, 250}, {0, 250}, {50, 50});
    hm.Add(met_vs_all_nu_py, "MET vs sum of all neutrinos in the event", {"MET py, [GeV]", "all nu py, [GeV]"}, {0, 250}, {0, 250}, {50, 50});
    hm.Add(jet_vs_quark_b, "Jet pt vs  quark pt: b", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_bbar, "Jet pt vs  quark pt: bbar", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_b_wnu, "Jet pt vs  quark pt: b (neutrino added)", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_bbar_wnu, "Jet pt vs  quark pt: bbar (neutrino added)", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_q1, "Jet pt vs  quark pt: q1", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_q2, "Jet pt vs  quark pt: q2", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_q1_wnu, "Jet pt vs  quark pt: q1 (neutrino added)", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(jet_vs_quark_q2_wnu, "Jet pt vs  quark pt: q2 (neutrino added)", {"Jet pt, [GeV]", "quark, [GeV]"}, {0, 350}, {0, 350}, {50, 50});
    hm.Add(hadW_vs_lepW, "Hadronic W mass vs Leptonic W mass", {"Hadronic W mass, [GeV]", "Leptonic W mass, [GeV]"}, {0, 120}, {0, 120}, {50, 50});

    std::cout << std::boolalpha;
    int X_mass = 0; 

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            ++has_bbWW_decay;
            X_mass = static_cast<int>(GenPart_mass[sig[SIG::X]]);

            // unsigned int n = std::count_if(GenJetAK4_pt, GenJetAK4_pt + nGenJetAK4, [](double pt){ return pt >= 10.0; });
            // assert(n == nGenJetAK4);

            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                ++has_only_emu;

                hm.Fill(hadW_vs_lepW, GenPart_mass[sig[SIG::HadWlast]], GenPart_mass[sig[SIG::LepWlast]]);

                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                GenJetData gen_jet_data{genjet_ak4, {GenJetAK4_partonFlavour, GenJetAK4_hadronFlavour}};

                // CUTS START
                if (!PassLeptonCut(l_p4))
                {
                    continue;
                }
                ++has_accept_lep;

                std::vector<int> selected_jets = PrimaryJetSelection(genjet_ak4);
                if (selected_jets.size() < 4)
                {
                    continue;
                }
                ++has_at_least_4_primary_jets;   

                // perform jet cleaning: reject jets overlapping with the lepton
                auto JetLepOverlap = [&l_p4, &genjet_ak4](int i)
                { 
                    TLorentzVector p = GetP4(genjet_ak4, i);
                    return p.DeltaR(l_p4) < DR_THRESH;
                };
                selected_jets.erase(std::remove_if(selected_jets.begin(), selected_jets.end(), JetLepOverlap), selected_jets.end());
                if (selected_jets.size() < 4)
                {
                    continue;
                }
                ++has_at_least_4_clean_jets;

                std::vector<TLorentzVector> b_jets;
                for (auto j: selected_jets)
                {
                    unsigned flav = static_cast<unsigned>(GenJetAK4_hadronFlavour[j]);
                    if (flav == 5)
                    {
                        b_jets.push_back(GetP4(genjet_ak4, j));
                    }
                }

                int num_bflav_jets = b_jets.size();
                if (num_bflav_jets < 1)
                {
                    continue;
                }
                ++has_at_least_1_bflav_jet;

                if (num_bflav_jets < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet;

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_B_GENJET_PT; };
                int num_bjets_passing_pt = std::count_if(b_jets.begin(), b_jets.end(), PtCut);
                if (num_bjets_passing_pt < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet_passing_pt;

                auto EtaCut = [](TLorentzVector const& p) { return std::abs(p.Eta()) < MAX_GENJET_ETA; };
                int num_bjets_passing_eta = std::count_if(b_jets.begin(), b_jets.end(), EtaCut);
                if (num_bjets_passing_eta < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet_passing_eta;

                int b_idx = sig[SIG::b];
                int bbar_idx = sig[SIG::bbar];
                int q1_idx = sig[SIG::q1];
                int q2_idx = sig[SIG::q2];

                // DO MATCHING HERE
                
                int b_match = Match(b_idx, genpart, genjet_ak4);
                int bbar_match = Match(bbar_idx, genpart, genjet_ak4);
                int q1_match = Match(q1_idx, genpart, genjet_ak4);
                int q2_match = Match(q2_idx, genpart, genjet_ak4);
        
                // skip event if matching failed
                std::vector<int> matches{b_match, bbar_match, q1_match, q2_match};
                std::vector<int> matches_copy = matches;
                // need to sort matches before doing unique
                // unique removes all except the first CONSECUTIVE duplicates 
                // (see https://en.cppreference.com/w/cpp/algorithm/unique)
                std::sort(matches.begin(), matches.end());
                bool negative_match = std::any_of(matches.begin(), matches.end(), [](int x) { return x == -1;});
                bool same_match = (std::unique(matches.begin(), matches.end()) != matches.end());
                if(negative_match || same_match)
                {
                    continue;
                } 
                else
                {
                    ++matched_events;

                    TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                    TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);
                    TLorentzVector bj_p4 = GetP4(genjet_ak4, b_match);
                    TLorentzVector bbarj_p4 = GetP4(genjet_ak4, bbar_match);

                    TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                    TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                    TLorentzVector lj1_p4 = GetP4(genjet_ak4, q1_match);
                    TLorentzVector lj2_p4 = GetP4(genjet_ak4, q2_match);

                    if (lj1_p4.DeltaR(l_p4) < DR_THRESH)
                    {
                        lj1_p4 -= l_p4;
                    }
                    if (lj2_p4.DeltaR(l_p4) < DR_THRESH)
                    {
                        lj2_p4 -= l_p4;
                    }

                    TLorentzVector genmet;
                    genmet.SetPtEtaPhiM(GenMET_pt, 0, GenMET_phi, 0);

                    if (!ConsistentMatch({lq1_p4, lj1_p4}, {lq2_p4, lj2_p4}))
                    {
                        ++inconsistent_light_jet_match;
                        continue;
                    }

                    if (!ConsistentMatch({bq_p4, bj_p4}, {bbarq_p4, bbarj_p4}))
                    {
                        ++inconsistent_b_jet_match;
                        continue;
                    }

                    ++accepted_events;                

                    #ifdef DEBUG
                        double leading_b_pt_ratio = (bq_p4.Pt() > bbarq_p4.Pt()) ? (bq_p4.Pt()/bj_p4.Pt()) : (bbarq_p4.Pt()/bbarj_p4.Pt());
                        double leading_l_pt_ratio = (lq1_p4.Pt() > lq2_p4.Pt()) ? (lq1_p4.Pt()/lj1_p4.Pt()) : (lq2_p4.Pt()/lj2_p4.Pt());

                        double subleading_b_pt_ratio = (bq_p4.Pt() > bbarq_p4.Pt()) ? (bbarq_p4.Pt()/bbarj_p4.Pt()) : (bq_p4.Pt()/bj_p4.Pt());
                        double subleading_l_pt_ratio = (lq1_p4.Pt() > lq2_p4.Pt()) ? (lq2_p4.Pt()/lj2_p4.Pt()) : (lq1_p4.Pt()/lj1_p4.Pt());

                        std::vector<double> ratios{ leading_b_pt_ratio, subleading_b_pt_ratio, leading_l_pt_ratio, subleading_l_pt_ratio };
                        if (std::any_of(ratios.begin(), ratios.end(), [](double x) { return x > MATCH_PT_THRESH; }))
                        {
                            MatchIndex mi = {{b_idx, b_match}, {bbar_idx, bbar_match}, {q1_idx, q1_match}, {q2_idx, q2_match}};
                            MatchKinematics mk = {genpart, genjet_ak4};
                            std::pair<int*, int*> ptrs = {GenPart_genPartIdxMother, GenPart_pdgId};
                            DrawEventMap(mk, mi, i, ptrs);
                        }
                    #endif

                    // determines if particle at index idx is a neutrino
                    auto NotNu = [&GenPart_pdgId](int idx) { return !IsNu(GenPart_pdgId[idx]); };

                    auto all_nus = GetFinalParticles(GenPart_genPartIdxMother, nGenPart);
                    all_nus.erase(std::remove_if(all_nus.begin(), all_nus.end(), NotNu), all_nus.end());

                    std::vector<int> b_neutrinos;
                    std::copy_if(all_nus.begin(), all_nus.end(), std::back_inserter(b_neutrinos), 
                                [&b_idx, &GenPart_genPartIdxMother](int idx) { return IsDescOf(idx, b_idx, GenPart_genPartIdxMother); });

                    std::vector<int> bbar_neutrinos;
                    std::copy_if(all_nus.begin(), all_nus.end(), std::back_inserter(bbar_neutrinos), 
                                [&bbar_idx, &GenPart_genPartIdxMother](int idx) { return IsDescOf(idx, bbar_idx, GenPart_genPartIdxMother); });

                    std::vector<int> q1_neutrinos;
                    std::copy_if(all_nus.begin(), all_nus.end(), std::back_inserter(q1_neutrinos), 
                                [&q1_idx, &GenPart_genPartIdxMother](int idx) { return IsDescOf(idx, q1_idx, GenPart_genPartIdxMother); });

                    std::vector<int> q2_neutrinos;
                    std::copy_if(all_nus.begin(), all_nus.end(), std::back_inserter(q2_neutrinos), 
                                [&q2_idx, &GenPart_genPartIdxMother](int idx) { return IsDescOf(idx, q2_idx, GenPart_genPartIdxMother); });

                    // 4-momentum generator
                    auto GenP4 = [&genpart](int idx) { return GetP4(genpart, idx); };

                    TLorentzVector tot_nu_p4 = std::transform_reduce(all_nus.begin(), all_nus.end(), TLorentzVector(), std::plus<TLorentzVector>(), GenP4);
                    TLorentzVector tot_nu_b_p4 = std::transform_reduce(b_neutrinos.begin(), b_neutrinos.end(), TLorentzVector(), std::plus<TLorentzVector>(), GenP4);
                    TLorentzVector tot_nu_bbar_p4 = std::transform_reduce(bbar_neutrinos.begin(), bbar_neutrinos.end(), TLorentzVector(), std::plus<TLorentzVector>(), GenP4);
                    TLorentzVector tot_nu_q1_p4 = std::transform_reduce(q1_neutrinos.begin(), q1_neutrinos.end(), TLorentzVector(), std::plus<TLorentzVector>(), GenP4);
                    TLorentzVector tot_nu_q2_p4 = std::transform_reduce(q2_neutrinos.begin(), q2_neutrinos.end(), TLorentzVector(), std::plus<TLorentzVector>(), GenP4);

                    hm.Fill(met_vs_all_nu, GenMET_pt, tot_nu_p4.Pt());
                    hm.Fill(met_vs_all_nu_px, genmet.Px(), tot_nu_p4.Px());
                    hm.Fill(met_vs_all_nu_py, genmet.Py(), tot_nu_p4.Py());

                    hm.Fill(jet_vs_quark_b, bj_p4.Pt(), bq_p4.Pt());
                    hm.Fill(jet_vs_quark_bbar, bbarj_p4.Pt(), bbarq_p4.Pt());
                    hm.Fill(higgs_mass_jets, (bj_p4 + bbarj_p4).M());
                    hm.Fill(w_mass_jets, (lj1_p4 + lj2_p4).M());
                    hm.Fill(pt_ratio_b, (bj_p4.Pt() - bq_p4.Pt())/bq_p4.Pt());
                    hm.Fill(pt_ratio_bbar, (bbarj_p4.Pt() - bbarq_p4.Pt())/bbarq_p4.Pt());
                    hm.Fill(pt_ratio_1, (lj1_p4.Pt() - lq1_p4.Pt())/lq1_p4.Pt());
                    hm.Fill(pt_ratio_2, (lj2_p4.Pt() - lq2_p4.Pt())/lq2_p4.Pt());

                    hm.Fill(jet_vs_quark_q1, lj1_p4.Pt(), lq1_p4.Pt());
                    hm.Fill(jet_vs_quark_q2, lj2_p4.Pt(), lq2_p4.Pt());

                    bj_p4 += tot_nu_b_p4;
                    bbarj_p4 += tot_nu_bbar_p4;
                    lj1_p4 += tot_nu_q1_p4;
                    lj2_p4 += tot_nu_q2_p4;

                    hm.Fill(jet_vs_quark_b_wnu, bj_p4.Pt(), bq_p4.Pt());
                    hm.Fill(jet_vs_quark_bbar_wnu, bbarj_p4.Pt(), bbarq_p4.Pt());
                    hm.Fill(jet_vs_quark_q1_wnu, lj1_p4.Pt(), lq1_p4.Pt());
                    hm.Fill(jet_vs_quark_q2_wnu, lj2_p4.Pt(), lq2_p4.Pt());
                    hm.Fill(higgs_mass_jets_wnu, (bj_p4 + bbarj_p4).M());
                    hm.Fill(w_mass_jets_wnu, (lj1_p4 + lj2_p4).M());
                    hm.Fill(pt_ratio_b_wnu, (bj_p4.Pt() - bq_p4.Pt())/bq_p4.Pt());
                    hm.Fill(pt_ratio_bbar_wnu, (bbarj_p4.Pt() - bbarq_p4.Pt())/bbarq_p4.Pt());
                    hm.Fill(pt_ratio_1_wnu, (lj1_p4.Pt() - lq1_p4.Pt())/lq1_p4.Pt());
                    hm.Fill(pt_ratio_2_wnu, (lj2_p4.Pt() - lq2_p4.Pt())/lq2_p4.Pt());
                }
            }
            else
            {
                ++has_tau;
            }
        }
    }

    std::cout << "Finished processing\n";
    hm.Draw();
    hm.DrawStack({higgs_mass_jets, higgs_mass_jets_wnu}, "Higgs mass: effect of adding neutrinos to b jets", "higgs_mass");
    hm.DrawStack({w_mass_jets, w_mass_jets_wnu}, "W mass: effect of adding neutrinos to light jets", "w_mass");
    hm.DrawStack({pt_ratio_b, pt_ratio_b_wnu}, "(Jet_pt - quark_pt)/quark_pt: effect of adding neutrinos to b jets", "b_jet_ratios");
    hm.DrawStack({pt_ratio_bbar, pt_ratio_bbar_wnu}, "(Jet_pt - quark_pt)/quark_pt: effect of adding neutrinos to bbar jets", "bbar_jet_ratios");
    hm.DrawStack({pt_ratio_1, pt_ratio_1_wnu}, "(Jet_pt - quark_pt)/quark_pt: effect of adding neutrinos to light jets", "light_1_ratios");
    hm.DrawStack({pt_ratio_2, pt_ratio_2_wnu}, "(Jet_pt - quark_pt)/quark_pt: effect of adding neutrinos to light jets", "light_2_ratios");
    hm.DrawStack({pt_ratio_1, pt_ratio_2}, "(Jet_pt - quark_pt)/quark_pt", "light_jet_ratios");
    hm.DrawStack({pt_ratio_1, pt_ratio_2, pt_ratio_b, pt_ratio_bbar}, "(Jet_pt - quark_pt)/quark_pt", "all_jet_ratios_nonu");
    hm.DrawStack({pt_ratio_1, pt_ratio_2, pt_ratio_b_wnu, pt_ratio_bbar_wnu}, "(Jet_pt - quark_pt)/quark_pt", "all_jet_ratios_wnu");

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::vector<double> rel_eff{ 100.0*has_bbWW_decay/nEvents,
                                 100.0*has_only_emu/has_bbWW_decay, 
                                 100.0*has_accept_lep/has_only_emu, 
                                 100.0*has_at_least_4_primary_jets/has_accept_lep, 
                                 100.0*has_at_least_4_clean_jets/has_at_least_4_primary_jets,
                                 100.0*has_at_least_1_bflav_jet/has_at_least_4_clean_jets, 
                                 100.0*has_at_least_2_bflav_jet/has_at_least_1_bflav_jet,
                                 100.0*has_at_least_2_bflav_jet_passing_pt/has_at_least_2_bflav_jet,
                                 100.0*has_at_least_2_bflav_jet_passing_eta/has_at_least_2_bflav_jet_passing_pt };
                    
    std::vector<double> abs_eff{ 100.0*has_bbWW_decay/nEvents,
                                 100.0*has_only_emu/nEvents, 
                                 100.0*has_accept_lep/nEvents, 
                                 100.0*has_at_least_4_primary_jets/nEvents, 
                                 100.0*has_at_least_4_clean_jets/nEvents,
                                 100.0*has_at_least_1_bflav_jet/nEvents,
                                 100.0*has_at_least_2_bflav_jet/nEvents, 
                                 100.0*has_at_least_2_bflav_jet_passing_pt/nEvents,
                                 100.0*has_at_least_2_bflav_jet_passing_eta/nEvents };

    std::vector<char const*> cut_labels{ "has bbWW", 
                                         "only e(mu)",
                                         "accepted lep",
                                         ">= 4 jets", 
                                         ">= 4 clean jets", 
                                         ">= 1 b jet",
                                         ">= 2 b jets",
                                         ">= 2 b jets pt > 20",
                                         ">= 2 b jets |eta| < 2.5" };

    int nx = cut_labels.size();
    auto rel_eff_hist = std::make_unique<TH1F>("rel_eff_hist", Form("Relative cut efficiency: %d GeV", X_mass), nx, 0, nx);
    auto abs_eff_hist = std::make_unique<TH1F>("abs_eff_hist", Form("Absolute cut efficiency: %d GeV", X_mass), nx, 0, nx);

    rel_eff_hist->SetStats(0);
    rel_eff_hist->SetFillStyle(3544);
    rel_eff_hist->SetFillColorAlpha(kBlue, 0.75);
    abs_eff_hist->SetStats(0);
    abs_eff_hist->SetFillStyle(3544);
    abs_eff_hist->SetFillColorAlpha(kBlue, 0.75);

    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetBottomMargin(0.15);

    for (int b = 1; b <= nx; ++b)
    {
        rel_eff_hist->SetBinContent(b, rel_eff[b-1]);
        rel_eff_hist->GetXaxis()->SetBinLabel(b, cut_labels[b-1]);
        abs_eff_hist->SetBinContent(b, abs_eff[b-1]);
        abs_eff_hist->GetXaxis()->SetBinLabel(b, cut_labels[b-1]);
    }

    rel_eff_hist->Draw();
    c1->SaveAs("rel_eff_hist.png");

    abs_eff_hist->Draw();
    c1->SaveAs("abs_eff_hist.png");

    std::cout << std::setprecision(3);
    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n"
              << "\tHave bbWW decay: " << has_bbWW_decay << "/" << nEvents << " (" << 100.0*has_bbWW_decay/nEvents << "%)\n"
              << "\tHave tau leptons: " << has_tau << "/" << nEvents << " (" << 100.0*has_tau/nEvents << "%)\n" 
              << "\tHave only electrons/muons: " << has_only_emu << "/" << nEvents << " (" << 100.0*has_only_emu/nEvents << "%)\n"
              << "\tHave lepton in acceptance region: " << has_accept_lep << "/" << has_only_emu << " (" << 100.0*has_accept_lep/has_only_emu << "%)\n"
              << "\tHave at least 4 primary jets: " << has_at_least_4_primary_jets << "/" << has_accept_lep << " (" << 100.0*has_at_least_4_primary_jets/has_accept_lep << "%)\n"
              << "\tHave at least 4 clean jets: " << has_at_least_4_clean_jets << "/" << has_at_least_4_primary_jets << " (" << 100.0*has_at_least_4_clean_jets/has_at_least_4_primary_jets << "%)\n"
              << "\tHave at least 1 b-flavored jet: " << has_at_least_1_bflav_jet << "/" << has_at_least_4_clean_jets << " (" << 100.0*has_at_least_1_bflav_jet/has_at_least_4_clean_jets << "%)\n"
              << "\tHave at least 2 b-flavored jets: " << has_at_least_2_bflav_jet << "/" << has_at_least_1_bflav_jet << " (" << 100.0*has_at_least_2_bflav_jet/has_at_least_1_bflav_jet << "%)\n"
              << "\tHave at least 2 b-flavored jets passing pt cut: " << has_at_least_2_bflav_jet_passing_pt << "/" << has_at_least_2_bflav_jet << " (" << 100.0*has_at_least_2_bflav_jet_passing_pt/has_at_least_2_bflav_jet << "%)\n"
              << "\tHave at least 2 b-flavored jets passing eta cut: " << has_at_least_2_bflav_jet_passing_eta << "/" << has_at_least_2_bflav_jet_passing_pt << " (" << 100.0*has_at_least_2_bflav_jet_passing_eta/has_at_least_2_bflav_jet_passing_pt << "%)\n";
    return 0;
}