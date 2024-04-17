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
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-450.root");
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

    // ak4 jets
    UInt_t          nGenJetAK4;
    Float_t         GenJetAK4_eta[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_mass[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_phi[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_pt[MAX_AK4_GENJET];   //[nGenJetAK4]
    Int_t           GenJetAK4_partonFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    UChar_t         GenJetAK4_hadronFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    
    // ak8 jets
    UInt_t          nGenJetAK8;
    Float_t         GenJetAK8_eta[MAX_AK8_GENJET];   //[nGenJetAK8]
    Float_t         GenJetAK8_mass[MAX_AK8_GENJET];   //[nGenJetAK8]
    Float_t         GenJetAK8_phi[MAX_AK8_GENJET];   //[nGenJetAK8]
    Float_t         GenJetAK8_pt[MAX_AK8_GENJET];   //[nGenJetAK8]
    Int_t           GenJetAK8_partonFlavour[MAX_AK8_GENJET];   //[nGenJetAK8]
    UChar_t         GenJetAK8_hadronFlavour[MAX_AK8_GENJET];   //[nGenJetAK8]

    // ak8 subjets
    UInt_t          nSubGenJetAK8;
    Float_t         SubGenJetAK8_eta[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    Float_t         SubGenJetAK8_mass[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    Float_t         SubGenJetAK8_phi[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]
    Float_t         SubGenJetAK8_pt[2*MAX_AK8_GENJET];   //[nSubGenJetAK8]

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
    myTree->SetBranchAddress("nGenJetAK8", &nGenJetAK8);
    myTree->SetBranchAddress("GenJetAK8_eta", &GenJetAK8_eta);
    myTree->SetBranchAddress("GenJetAK8_mass", &GenJetAK8_mass);
    myTree->SetBranchAddress("GenJetAK8_phi", &GenJetAK8_phi);
    myTree->SetBranchAddress("GenJetAK8_pt", &GenJetAK8_pt);
    myTree->SetBranchAddress("GenJetAK8_partonFlavour", &GenJetAK8_partonFlavour);
    myTree->SetBranchAddress("GenJetAK8_hadronFlavour", &GenJetAK8_hadronFlavour);

    // ak8 hen subjets
    myTree->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8);
    myTree->SetBranchAddress("SubGenJetAK8_eta", &SubGenJetAK8_eta);
    myTree->SetBranchAddress("SubGenJetAK8_mass", &SubGenJetAK8_mass);
    myTree->SetBranchAddress("SubGenJetAK8_phi", &SubGenJetAK8_phi);
    myTree->SetBranchAddress("SubGenJetAK8_pt", &SubGenJetAK8_pt);

    myTree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    myTree->SetBranchAddress("GenMET_pt", &GenMET_pt);

    int nEvents = myTree->GetEntries();

    // counters

    int valid_sig = 0;
    int have_matchable_jets = 0;
    int matched_events = 0;
    int failed_genjet_cut = 0;
    int failed_lepton_cut = 0;
    int not_isolated_lepton = 0;
    int inconsistent_light_jet_match = 0;
    int inconsistent_b_jet_match = 0;
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

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        // std::cout << "Event " << i << ":\n";
        // std::cout << "\tAK8 jets pt: ";
        // for (int i = 0; i < static_cast<int>(nGenJetAK8); ++i)
        // {
        //     std::cout << GenJetAK8_pt[i] << " ";
        // }
        // std::cout << "\n";
        // std::cout << "\tAK8 subjetsjets pt: ";
        // for (int i = 0; i < static_cast<int>(nSubGenJetAK8); ++i)
        // {
        //     std::cout << GenJetAK8_pt[i] << " ";
        // }
        
        // if (std::cin.get())
        // {
        //     std::cout << "\n";
        //     continue;
        // }
        // std::cout << "===================================\n";

        if (!sig.empty())
        {
            if (CheckSignal(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                ++valid_sig;

                hm.Fill(hadW_vs_lepW, GenPart_mass[sig[SIG::HadWlast]], GenPart_mass[sig[SIG::LepWlast]]);

                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};
                KinematicData genjet_ak8{GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, GenJetAK8_mass, static_cast<int>(nGenJetAK8)};
                KinematicData subgenjet_ak8{SubGenJetAK8_pt, SubGenJetAK8_eta, SubGenJetAK8_phi, SubGenJetAK8_mass, static_cast<int>(nSubGenJetAK8)};

                auto matchable_jets = GetMatchableJets(genjet_ak4);
                int n_matchable_jets = matchable_jets.size();
                if (n_matchable_jets >= 4)
                {
                    ++have_matchable_jets;
                }
                // else
                // {
                //     continue;
                // }

                int b_idx = sig[SIG::b];
                int bbar_idx = sig[SIG::bbar];
                int q1_idx = sig[SIG::q1];
                int q2_idx = sig[SIG::q2];

                TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);
                TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                
                // std::cout << "Event " << i << ":\n";
                // std::cout << "matching " << GenPart_pdgId[b_idx] << "(pt = " << GenPart_pt[b_idx] << "):\n";
                int b_match = Match(b_idx, genpart, genjet_ak4);
                // std::cout << "matching " << GenPart_pdgId[bbar_idx] << "(pt = " << GenPart_pt[bbar_idx] << "):\n";
                int bbar_match = Match(bbar_idx, genpart, genjet_ak4);
                // std::cout << "matching " << GenPart_pdgId[q1_idx] << "(pt = " << GenPart_pt[q1_idx] << "):\n";
                int q1_match = Match(q1_idx, genpart, genjet_ak4);
                // std::cout << "matching " << GenPart_pdgId[q2_idx] << "(pt = " << GenPart_pt[q2_idx] << "):\n";
                int q2_match = Match(q2_idx, genpart, genjet_ak4);
                // std::cout << "======================================\n";

                // if (std::cin.get())
                // {
                //     continue;
                // }

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
                    // std::cout << "Event " << i << ":\n";
                    // std::cout << "\tmatches: ";
                    // for (auto m: matches_copy)
                    // {
                    //     std::cout << m << " ";
                    // }
                    // std::cout << "\n";
                    // std::cout << "\tmin_dR = " << MinDeltaR({bq_p4, bbarq_p4, lq1_p4, lq2_p4}) << "\n";
                    // std::cout << "\tdR(q, q) = " << lq1_p4.DeltaR(lq2_p4) << "\n";
                    // std::cout << "\tq1_pt = " << lq1_p4.Pt() << "\n";
                    // std::cout << "\tq2_pt = " << lq2_p4.Pt() << "\n";
                    // std::cout << "\tHadWlast_pt = " << GenPart_pt[sig[SIG::HadWlast]] << "\n\n";

                    // std::cout << "\tdR(b, b) = " << bq_p4.DeltaR(bbarq_p4) << "\n";
                    // std::cout << "\tb_pt = " << bq_p4.Pt() << "\n";
                    // std::cout << "\tbbar_pt = " << bbarq_p4.Pt() << "\n";
                    // std::cout << "\tH1_pt = " << GenPart_pt[sig[SIG::H1]] << "\n";
                    // std::cout << "\tH2_pt = " << GenPart_pt[sig[SIG::H2]] << "\n\n";

                    // std::cout << "\tnGenJetAK8 = " << nGenJetAK8 << "\n";
                    // if (nGenJetAK8 > 0)
                    // {
                    //     std::cout << "\tGenJetAK8_pt: ";
                    //     for (int j = 0; j < static_cast<int>(nGenJetAK8); ++j)
                    //     {
                    //         std::cout << GenJetAK8_pt[j] << " ";
                    //     }
                    //     std::cout << "\n";
                    // }

                    // std::cout << "==========================\n";
                    continue;

                    // if (std::cin.get())
                    // {
                    //     continue;
                    // }
                } 
                else
                {
                    ++matched_events;

                    // TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                    // TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);
                    TLorentzVector bj_p4 = GetP4(genjet_ak4, b_match);
                    TLorentzVector bbarj_p4 = GetP4(genjet_ak4, bbar_match);

                    // TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                    // TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                    TLorentzVector lj1_p4 = GetP4(genjet_ak4, q1_match);
                    TLorentzVector lj2_p4 = GetP4(genjet_ak4, q2_match);

                    TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                    TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                    TLorentzVector genmet;
                    genmet.SetPtEtaPhiM(GenMET_pt, 0, GenMET_phi, 0);

                    // apply cuts here
                    std::vector<TLorentzVector> jets{lj1_p4, lj2_p4, bj_p4, bbarj_p4};

                    if (!PassGenJetCut(jets))
                    {
                        ++failed_genjet_cut;
                        continue;
                    }  

                    if (!PassLeptonCut(l_p4))
                    {
                        ++failed_lepton_cut;
                        continue;
                    } 

                    if (!IsIsolatedLepton(l_p4, jets))
                    {
                        ++not_isolated_lepton;
                        continue;
                    }

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

                    // if (lj1_p4.Pt()/lq1_p4.Pt() > 5)
                    // {
                    //     std::cout << "Event " << i << ": light jet 1 too large\n";
                    //     std::cout << "q1 = "; Print(lq1_p4);
                    //     std::cout << "q2 = "; Print(lq2_p4);
                    //     std::cout << "b = "; Print(bq_p4);
                    //     std::cout << "bbar = "; Print(bbarq_p4);
                    //     std::cout << "================================\n";
                    // }
                    // if (lj2_p4.Pt()/lq2_p4.Pt() > 5)
                    // {
                    //     std::cout << "Event " << i << ": light jet 2 too large\n";
                    //     std::cout << "q1 = "; Print(lq1_p4);
                    //     std::cout << "q2 = "; Print(lq2_p4);
                    //     std::cout << "b = "; Print(bq_p4);
                    //     std::cout << "bbar = "; Print(bbarq_p4);
                    //     std::cout << "================================\n";
                    // }
                }
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

    std::cout << std::setprecision(3);
    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n" 
              << "\tAre signal events: " << valid_sig << "/" << nEvents << " (" << 100.0*valid_sig/nEvents << "%)\n"
              << "\tHave matchable jets: " << have_matchable_jets << "/" << valid_sig << " (" << 100.0*have_matchable_jets/valid_sig << "%)\n"
              << "\tSuccessfully matched all jets: " << matched_events << "/" << have_matchable_jets << " (" << 100.0*matched_events/have_matchable_jets << "%)\n"
              << "\tFailed genjet_ak4 cut: " << failed_genjet_cut << "/" << matched_events << " (" << 100.0*failed_genjet_cut/matched_events << "%)\n"
              << "\tFailed lepton cut: " << failed_lepton_cut << "/" << matched_events << " (" << 100.0*failed_lepton_cut/matched_events << "%)\n"
              << "\tNot isolated lepton: " << not_isolated_lepton << "/" << matched_events << " (" << 100.0*not_isolated_lepton/matched_events << "%)\n"
              << "\tInconsistet light jet matching: " << inconsistent_light_jet_match << "/" << matched_events << " (" << 100.0*inconsistent_light_jet_match/matched_events << "%)\n"
              << "\tInconsistet b jet matching: " << inconsistent_b_jet_match << "/" << matched_events << " (" << 100.0*inconsistent_b_jet_match/matched_events << "%)\n"
              << "\tEvents passed to HME = " << accepted_events << "\n";
    return 0;
}