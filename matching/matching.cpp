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
#include "TRandom3.h"

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
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-450.root");
    // TFile *myFile = TFile::Open("../../Di-Higgs/code/matching/GluGluToRadionToHHTo2B2WToLNu2J_M-1000.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    auto file_pdf = std::make_unique<TFile>("pdfs_hadW.root", "READ");
    auto pdf_hadW_offshell = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_hadW_offshell")));
    auto pdf_hadW_onshell = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_hadW_onshell")));

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

    // ak4 gen jets
    myTree->SetBranchAddress("nGenJet", &nGenJetAK4);
    myTree->SetBranchAddress("GenJet_eta", &GenJetAK4_eta);
    myTree->SetBranchAddress("GenJet_mass", &GenJetAK4_mass);
    myTree->SetBranchAddress("GenJet_phi", &GenJetAK4_phi);
    myTree->SetBranchAddress("GenJet_pt", &GenJetAK4_pt);
    myTree->SetBranchAddress("GenJet_partonFlavour", &GenJetAK4_partonFlavour);
    myTree->SetBranchAddress("GenJet_hadronFlavour", &GenJetAK4_hadronFlavour);

    myTree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    myTree->SetBranchAddress("GenMET_pt", &GenMET_pt);

    int nEvents = myTree->GetEntries();
    
    TRandom3 rg;
    rg.SetSeed(42);

    int WtoB = 0;
    int tot = 0;
    int selected_events = 0;
    int both_have_match = 0;
    int same_best_match = 0;
    int only_q1_match = 0;
    int only_q2_match = 0;

    int dijet_match_fail = 0;

    // hists
    HistManager hm;

    std::string hadW_mass_strategy_1("hadW_mass_strategy_1");
    hm.Add(hadW_mass_strategy_1, "Mass of HadW from jets chosen by minimizing metric", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_matched_dijet("hadW_mass_matched_dijet");
    hm.Add(hadW_mass_matched_dijet, "Mass of HadW from matching dijets to HadW", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_gen("hadW_mass_gen");
    hm.Add(hadW_mass_gen, "Gen mass of HadW", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_best_offshell_mass("hadW_best_offshell_mass");
    hm.Add(hadW_best_offshell_mass, "Best offshell HadW mass", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_best_onshell_mass("hadW_best_onshell_mass");
    hm.Add(hadW_best_onshell_mass, "Best onshell HadW mass", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_strategy_1_all("hadW_mass_strategy_1_all");
    hm.Add(hadW_mass_strategy_1_all, "Mass of HadW from jets chosen by minimizing metric for all accepted events", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_strategy_2_all("hadW_mass_strategy_2_all");
    hm.Add(hadW_mass_strategy_2_all, "Mass of HadW when choosing best onshell and offshell pairs separately", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_strategy_2("hadW_mass_strategy_2");
    hm.Add(hadW_mass_strategy_2, "Mass of HadW from jets chosen by selecting best offshel and onshell pairs separately and choosing best", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_mass_matched_single("hadW_mass_matched_single");
    hm.Add(hadW_mass_matched_single, "Mass of HadW from jets chosen by mathcing jets to quarks", {"HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_offshell_matched("hadW_offshell_matched");
    hm.Add(hadW_offshell_matched, "Visible mass of offshell HadW", {"offshell HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_onshell_matched("hadW_onshell_matched");
    hm.Add(hadW_onshell_matched, "Visible mass of onshell HadW", {"onshell HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_offshell_gen("hadW_offshell_gen");
    hm.Add(hadW_offshell_gen, "Gen mass of offshell HadW", {"offshell HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    std::string hadW_onshell_gen("hadW_onshell_gen");
    hm.Add(hadW_onshell_gen, "Gen mass of onshell HadW", {"onshell HadW mass, [GeV]", "Count"}, {0, 200}, 100);

    // std::string pdf_hadW_offshell("pdf_hadW_offshell");
    // hm.Add(pdf_hadW_offshell, "PDF for visible mass of offshell HadW", {"offshell HadW mass, [GeV]", "Count"}, {0, 300}, 300);

    // std::string pdf_hadW_onshell("pdf_hadW_onshell");
    // hm.Add(pdf_hadW_onshell, "PDF for visible mass of onshell HadW", {"onshell HadW mass, [GeV]", "Count"}, {0, 300}, 300);

    // std::string hadW_mass_choose("hadW_mass_choose");
    // hm.Add(hadW_mass_choose, "Mass of HadW from jets chosen by minimizing angle between them and |m(jj) - mW|", {"HadW mass, [GeV]", "Count"}, {0, 300}, 100);

    std::string number_of_matches("number_of_matches");
    hm.Add(number_of_matches, "Number of jets matched  to q1 and q2", {"N", "Count"}, {-0.5, 6.5}, 7);

    // std::string q1_pt_vs_q2_pt_q1_match("q1_pt_vs_q2_pt_q1_match");
    // hm.Add(q1_pt_vs_q2_pt_q1_match, "q1 pt vs q2 pt when q1 has match", {"q1 pt, [GeV]", "q2 pt, [GeV]"}, {0, 200}, {0, 200}, {100, 100});

    // std::string q1_pt_vs_q2_pt_q2_match("q1_pt_vs_q2_pt_q2_match");
    // hm.Add(q1_pt_vs_q2_pt_q2_match, "q1 pt vs q2 pt when q2 has match", {"q1 pt, [GeV]", "q2 pt, [GeV]"}, {0, 200}, {0, 200}, {100, 100});

    std::string pt_of_quark_wo_match("pt_of_quark_wo_match");
    hm.Add(pt_of_quark_wo_match, "Pt of quark which doesn't have match within dR < 0.4", {"Pt, [GeV]", "Count"}, {0, 200}, 100);

    std::string min_dr_of_quark_wo_match("min_dr_of_quark_wo_match");
    hm.Add(min_dr_of_quark_wo_match, "Min dR between quark without match and selected jets", {"dR", "Count"}, {0, 6}, 100);

    std::string best_offshell_mass_vs_best_onshell_mass("best_offshell_mass_vs_best_onshell_mass");
    hm.Add(best_offshell_mass_vs_best_onshell_mass, "Best mass of offshell pair vs best mass of onshell pair", {"Offshell mass, [GeV]", "Onshell mass, [GeV]"}, {0, 200}, {0, 200}, {100, 100});

    std::string number_of_light_jets("number_of_light_jets");
    hm.Add(number_of_light_jets, "Number of selected light jets", {"N", "Count"}, {-0.5, 10.5}, 11);

    std::cout << std::boolalpha;

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                ++tot;
                if (std::abs(GenPart_pdgId[sig[SIG::q1]]) == 5 || std::abs(GenPart_pdgId[sig[SIG::q2]]) == 5)
                {
                    ++WtoB;
                }

                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                TLorentzVector bq1_p4 = GetP4(genpart, sig[SIG::b]);
                TLorentzVector bq2_p4 = GetP4(genpart, sig[SIG::bbar]);

                TLorentzVector q1_p4 = GetP4(genpart, sig[SIG::q1]);
                TLorentzVector q2_p4 = GetP4(genpart, sig[SIG::q2]);

                if (q1_p4.DeltaR(q2_p4) < DR_THRESH)
                {
                    continue;
                }

                if (q1_p4.Pt() < MIN_GENJET_PT || q2_p4.Pt() < MIN_GENJET_PT)
                {
                    continue;
                }

                if (std::abs(q1_p4.Eta()) > MAX_GENJET_ETA || std::abs(q2_p4.Eta()) > MAX_GENJET_ETA)
                {
                    continue;
                }

                if (q1_p4.DeltaR(l_p4) < DR_THRESH || q2_p4.DeltaR(l_p4) < DR_THRESH)
                {
                    continue;
                }

                TLorentzVector Hbb_p4 = GetP4(genpart, sig[SIG::H_bb]);
            
                // CUTS START
                if (!PassLeptonCut(l_p4))
                {
                    continue;
                }

                std::vector<int> selected_jets = PrimaryJetSelection(genjet_ak4);
                if (selected_jets.size() < 4)
                {
                    continue;
                }

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

                std::vector<TLorentzVector> b_jets, light_jets;
                std::vector<int> light_jet_indices;
                for (auto j: selected_jets)
                {
                    unsigned flav = static_cast<unsigned>(GenJetAK4_hadronFlavour[j]);
                    if (flav == 5)
                    {
                        b_jets.push_back(GetP4(genjet_ak4, j));
                    }
                    else
                    {
                        light_jets.push_back(GetP4(genjet_ak4, j));
                        light_jet_indices.push_back(j);
                    }
                }

                int num_bflav_jets = b_jets.size();
                if (num_bflav_jets < 1)
                {
                    continue;
                }

                if (num_bflav_jets < 2)
                {
                    continue;
                }

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_GENJET_PT; };
                int num_bjets_passing_pt = std::count_if(b_jets.begin(), b_jets.end(), PtCut);
                if (num_bjets_passing_pt < 2)
                {
                    continue;
                }

                auto EtaCut = [](TLorentzVector const& p) { return std::abs(p.Eta()) < MAX_GENJET_ETA; };
                int num_bjets_passing_eta = std::count_if(b_jets.begin(), b_jets.end(), EtaCut);
                if (num_bjets_passing_eta < 2)
                {
                    continue;
                }

                if (light_jets.size() < 2)
                {
                    continue;
                }

                TLorentzVector bj1_p4, bj2_p4;
                if (b_jets.size() == 2)
                {
                    bj1_p4 = b_jets[0].DeltaR(bq1_p4) < b_jets[1].DeltaR(bq1_p4) ? b_jets[0] : b_jets[1];
                    bj2_p4 = b_jets[0].DeltaR(bq2_p4) < b_jets[1].DeltaR(bq2_p4) ? b_jets[0] : b_jets[1];
                }
                else
                {
                    continue;
                }

                if (bj1_p4 == bj2_p4)
                {
                    continue;
                }

                ++selected_events;

                hm.Fill(number_of_light_jets, light_jets.size());

                auto [mean_ang, std_ang] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return p1.Angle(p2.Vect()); });
                // auto [mean_dr, std_dr] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return p1.DeltaR(p2); });
                auto [mean_dm, std_dm] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return 80.0 - (p1 + p2).M(); });
                // auto [mean_pt, std_pt] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return (p1 + p2).Pt(); });

                // auto metric = [&mean_dr, &std_dr, &mean_dm, &std_dm](TLorentzVector const& p1, TLorentzVector const& p2)
                auto metric = [&mean_ang, &std_ang, &mean_dm, &std_dm](TLorentzVector const& p1, TLorentzVector const& p2)
                {
                    double ang = p1.Angle(p2.Vect()); 
                    double z_ang = (ang - mean_ang)/std_ang;
                    
                    // double dr = p1.DeltaR(p2);
                    // double z_dr = (dr - mean_dr)/std_dr;

                    double dm = 80.0 - (p1 + p2).M();
                    double z_dm = (dm - mean_dm)/std_dm;

                    // double pt = (p1 + p2).Pt();
                    // double z_pt = (pt - mean_pt)/std_pt;

                    return z_ang*z_ang + z_dm*z_dm;
                    // return z_dr*z_dr + z_dm*z_dm;
                };

                std::vector<std::pair<int, int>> pairs;
                for (size_t i = 0; i < light_jets.size(); ++i)
                {
                    for (size_t j = i + 1; j < light_jets.size(); ++j)
                    {
                        pairs.emplace_back(i, j);
                    }
                }

                auto JetComparator = [&metric, &light_jets](std::pair<int, int> const& p1, std::pair<int, int> const& p2)
                {
                    return metric(light_jets[p1.first], light_jets[p1.second]) < metric(light_jets[p2.first], light_jets[p2.second]);
                };
                std::sort(pairs.begin(), pairs.end(), JetComparator);

                std::pair<int, int> bp{-1, -1};
                for (auto const& p: pairs)
                {
                    double mass = (l_p4 + light_jets[p.first] + light_jets[p.second]).M();
                    if (mass < 125.0)
                    {
                        bp = p;
                        break;
                    }
                }

                if (bp.first == -1)
                {
                    bp = pairs[0];
                }
                
                // auto [bi1, bi2] = ChooseBestPair(light_jets, metric);
                auto [bi1, bi2] = bp;
                double bp_mass = (light_jets[bi1] + light_jets[bi2]).M();
                hm.Fill(hadW_mass_strategy_1_all, bp_mass);

                auto best_onshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(80.0 - (p1 + p2).M()); });
                auto best_offshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(40.0 - (p1 + p2).M()); });

                double best_onshell_mass = (light_jets[best_onshell_pair.first] + light_jets[best_onshell_pair.second]).M();
                double best_offshell_mass = (light_jets[best_offshell_pair.first] + light_jets[best_offshell_pair.second]).M();
                hm.Fill(hadW_mass_strategy_2_all, rg.Uniform(0, 1) > 0.27 ? best_onshell_mass : best_offshell_mass);
                // hm.FillWeighted(hadW_mass_strategy_2_all, best_offshell_mass, 0.5);
                // hm.FillWeighted(hadW_mass_strategy_2_all, best_onshell_mass, 0.5);

                if (light_jets.size() > 2)
                {
                    hm.Fill(best_offshell_mass_vs_best_onshell_mass, best_offshell_mass, best_onshell_mass);
                }

                auto mq1 = Matches(sig[SIG::q1], genpart, genjet_ak4);
                auto mq2 = Matches(sig[SIG::q2], genpart, genjet_ak4);
                hm.Fill(number_of_matches, mq1.size() + mq2.size());

                // int best_q1_match = Match(sig[SIG::q1], genpart, genjet_ak4);
                // int best_q2_match = Match(sig[SIG::q2], genpart, genjet_ak4);
                int best_q1_match = Match(q1_p4, light_jets);
                int best_q2_match = Match(q2_p4, light_jets);

                TLorentzVector hadW = GetP4(genpart, sig[SIG::HadWlast]);
                TLorentzVector matched_dijet = MatchDijet(hadW, light_jets);
                hm.Fill(hadW_mass_matched_dijet, matched_dijet.M());
                // if (matched_dijet.Pt() > 0.0)
                // {
                //     hm.Fill(hadW_mass_matched_dijet, matched_dijet.M());
                // }
                // else 
                // {
                //     ++dijet_match_fail;
                // }

                if (best_q1_match != -1 && best_q2_match == -1)
                {
                    ++only_q1_match;
                    // hm.Fill(q1_pt_vs_q2_pt_q1_match, q1_p4.Pt(), q2_p4.Pt());
                    hm.Fill(pt_of_quark_wo_match, q2_p4.Pt());
                    hm.Fill(min_dr_of_quark_wo_match, MinDeltaR(q2_p4, light_jets));

                    // find second match by mass
                    int second_match = FindSecondMatch(best_q1_match, hadW, light_jets);
                    hm.Fill(hadW_mass_matched_single, (light_jets[best_q1_match] + light_jets[second_match]).M());
                }

                if (best_q1_match == -1 && best_q2_match != -1)
                {
                    ++only_q2_match;
                    // hm.Fill(q1_pt_vs_q2_pt_q2_match, q1_p4.Pt(), q2_p4.Pt());
                    hm.Fill(pt_of_quark_wo_match, q1_p4.Pt());
                    hm.Fill(min_dr_of_quark_wo_match, MinDeltaR(q1_p4, light_jets));

                    int second_match = FindSecondMatch(best_q2_match, hadW, light_jets);
                    hm.Fill(hadW_mass_matched_single, (light_jets[best_q2_match] + light_jets[second_match]).M());
                }

                if (best_q1_match != -1 && best_q2_match != -1)
                {
                    ++both_have_match;
                    if (best_q1_match == best_q2_match)
                    {
                        ++same_best_match;
                        hm.Fill(hadW_mass_matched_single, light_jets[best_q1_match].M());
                    }
                    else 
                    {
                        // TLorentzVector j1_p4 = GetP4(genjet_ak4, best_q1_match);
                        // TLorentzVector j2_p4 = GetP4(genjet_ak4, best_q2_match);
                        TLorentzVector j1_p4 = light_jets[best_q1_match];
                        TLorentzVector j2_p4 = light_jets[best_q2_match];
                        double hadW_mass = (j1_p4 + j2_p4).M();
                        hm.Fill(hadW_mass_matched_single, hadW_mass);

                        double hadW_gen = GenPart_mass[sig[SIG::HadWlast]];
                        double lepW_gen = GenPart_mass[sig[SIG::LepWlast]];
                        if (hadW_gen > lepW_gen)
                        {
                            hm.Fill(hadW_onshell_matched, hadW_mass);
                            hm.Fill(hadW_onshell_gen, hadW_gen);
                            // hm.Fill(pdf_hadW_onshell, hadW_mass);
                        }
                        else
                        {
                            hm.Fill(hadW_offshell_matched, hadW_mass);
                            hm.Fill(hadW_offshell_gen, hadW_gen);
                            // hm.Fill(pdf_hadW_offshell, hadW_mass);
                        }
                        hm.Fill(hadW_mass_gen, hadW_gen);

                        auto [mean_ang, std_ang] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return p1.Angle(p2.Vect()); });
                        auto [mean_dm, std_dm] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return 80.0 - (p1 + p2).M(); });
                        // auto [mean_pt, std_pt] = CalcJetPairStats(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return (p1 + p2).Pt(); });

                        auto metric = [&mean_ang, &std_ang, &mean_dm, &std_dm](TLorentzVector const& p1, TLorentzVector const& p2)
                        {
                            double ang = p1.Angle(p2.Vect()); 
                            double z_ang = (ang - mean_ang)/std_ang;

                            double dm = 80.0 - (p1 + p2).M();
                            double z_dm = (dm - mean_dm)/std_dm;

                            // double pt = (p1 + p2).Pt();
                            // double z_pt = (pt - mean_pt)/std_pt;

                            return z_ang*z_ang + z_dm*z_dm;
                        };

                        std::vector<std::pair<int, int>> pairs;
                        for (size_t i = 0; i < light_jets.size(); ++i)
                        {
                            for (size_t j = i + 1; j < light_jets.size(); ++j)
                            {
                                pairs.emplace_back(i, j);
                            }
                        }

                        auto JetComparator = [&metric, &light_jets](std::pair<int, int> const& p1, std::pair<int, int> const& p2)
                        {
                            return metric(light_jets[p1.first], light_jets[p1.second]) < metric(light_jets[p2.first], light_jets[p2.second]);
                        };
                        std::sort(pairs.begin(), pairs.end(), JetComparator);

                        // std::pair<int, int> bp{-1, -1};
                        bp = std::pair<int, int>(-1, -1);
                        for (auto const& p: pairs)
                        {
                            double mass = (l_p4 + light_jets[p.first] + light_jets[p.second]).M();
                            if (mass < 125.0)
                            {
                                bp = p;
                                break;
                            }
                        }

                        if (bp.first == -1)
                        {
                            bp = pairs[0];
                        }

                        auto [bi1, bi2] = bp;
                        double bp_mass = (light_jets[bi1] + light_jets[bi2]).M();
                        hm.Fill(hadW_mass_strategy_1, bp_mass);

                        best_onshell_mass = (light_jets[best_onshell_pair.first] + light_jets[best_onshell_pair.second]).M();
                        // double onshell_weight = pdf_hadW_onshell->GetBinContent(pdf_hadW_onshell->FindBin(best_onshell_mass));
                        hm.Fill(hadW_best_onshell_mass, best_onshell_mass);

                        best_offshell_mass = (light_jets[best_offshell_pair.first] + light_jets[best_offshell_pair.second]).M();
                        // double offshell_weight = pdf_hadW_offshell->GetBinContent(pdf_hadW_offshell->FindBin(best_offshell_mass));
                        hm.Fill(hadW_best_offshell_mass, best_offshell_mass);

                        // hm.Fill(best_offshell_mass_vs_best_onshell_mass, best_offshell_mass, best_onshell_mass);

                        // double best_mass = onshell_weight > offshell_weight ? best_onshell_mass : best_offshell_mass;
                        double best_mass = rg.Uniform(0, 1) > 0.27 ? best_onshell_mass : best_offshell_mass;
                        hm.Fill(hadW_mass_strategy_2, best_mass);
                        // hm.FillWeighted(hadW_mass_strategy_2, best_offshell_mass, 0.5);
                        // hm.FillWeighted(hadW_mass_strategy_2, best_onshell_mass, 0.5);
                    }
                }
            }    
        }
    }

    hm.Draw();

    hm.DrawStack({hadW_mass_matched_single, hadW_mass_strategy_1, hadW_mass_strategy_2}, "Comparison of different strategies to choose jets for W->qq with matching one jet to one quark", "strat_choice_vs_single_matching_cmp");
    hm.DrawStack({hadW_mass_matched_dijet, hadW_mass_strategy_1, hadW_mass_strategy_2}, "Comparison of different strategies to choose jets for W->qq with pair of jets to HadW", "strat_choice_vs_pair_matching_cmp");
    hm.DrawStack({hadW_mass_matched_dijet, hadW_mass_matched_single, hadW_mass_gen}, "Gen HadW mass vs different strategies to match jets", "matching_cmp");
    hm.DrawStack({hadW_mass_gen, hadW_mass_strategy_1, hadW_mass_strategy_2}, "Gen HadW mass vs different strategies to choose jets", "gen_vs_strat_cmp");
    // hm.DrawStack({hadW_best_onshell_mass, hadW_onshell_gen}, "Best onshell mass vs Gen onshell mass", "onshell_mass_cmp");
    // hm.DrawStack({hadW_best_offshell_mass, hadW_offshell_gen}, "Best offshell mass vs Gen offshell mass", "offshell_mass_cmp");

    hm.DrawStack({hadW_mass_matched_dijet, hadW_mass_matched_single}, "Comparison of different matching strategies", "matching_strat_cmp");
    hm.DrawStack({hadW_mass_strategy_1_all, hadW_mass_strategy_2_all}, "Comparison of different strategies to choose jets for all selected events", "choosing_strat_all_cmp");

    std::cout << "W decays to b: " << WtoB << "\n";
    std::cout << "Signal events: " << tot << "\n";

    std::cout << "Selected signal events: " << selected_events << "\n";
    std::cout << "Only q1 has match: " << only_q1_match << "\n";
    std::cout << "Only q2 has match: " << only_q2_match << "\n";
    std::cout << "Both light quarks have match: " << both_have_match << "\n";
    std::cout << "Both light quarks have same match: " << same_best_match << "\n";
    std::cout << "Dijet matching failed: " << dijet_match_fail << "\n";

    // hm.Write1D("pdfs_hadW.root", {pdf_hadW_offshell, pdf_hadW_onshell});

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing, time elapsed " << elapsed.count() << " s\n";

    return 0;
}