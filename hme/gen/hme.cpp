#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "MatchingTools.hpp"
#include "HistManager.hpp"
#include "EstimatorTools.hpp"

static constexpr int MAX_AK4_GENJET = 35;
static constexpr int MAX_GENPART = 350;

double InterquantileRange(std::unique_ptr<TH1F> const& h)
{
    int const nq = 50;
    double xq[nq];  // position where to compute the quantiles in [0,1]
    double yq[nq];  // array to contain the quantiles
    for (int i = 0; i < nq; ++i) 
    {
        xq[i] = static_cast<double>(i+1)/nq;
    }
    h->GetQuantiles(nq, yq, xq);
    return yq[41] - yq[7];
}

int main()
{
    auto start = std::chrono::system_clock::now();
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-450.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    auto file_pdf = std::make_unique<TFile>("pdfs.root", "READ");
    auto lead_bjet_pt_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("lead_bjet_pt_pdf")));
    auto pdf2d = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get("bjet_2d_pt_pdf")));

    TRandom3 rg;
    rg.SetSeed(0);

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

    // hists
    HistManager hm;

    std::string hme_mass("hme_mass");
    hm.Add(hme_mass, "HME X->HH mass with ideal H->WW", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string hme_mass_ideal_nu("hme_mass_ideal_nu");
    hm.Add(hme_mass_ideal_nu, "HME X->HH mass with ideal nu", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string hme_mass_pdf2d_ideal_nu("hme_mass_pdf2d_ideal_nu");
    hm.Add(hme_mass_pdf2d_ideal_nu, "HME X->HH mass with ideal nu and 2d PDF", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string hme_ideal_hadW("hme_ideal_hadW");
    hm.Add(hme_ideal_hadW, "HME X->HH mass with ideal W->qq", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string mass_ideal_HWW("mass_ideal_HWW");
    hm.Add(mass_ideal_HWW, "X->HH mass with ideal H->WW", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::cout << std::boolalpha;
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                int b_idx = sig[SIG::b];
                int bbar_idx = sig[SIG::bbar];
                int q1_idx = sig[SIG::q1];
                int q2_idx = sig[SIG::q2];

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
                if (num_bflav_jets < 2)
                {
                    continue;
                }

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_B_GENJET_PT; };
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

                TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);
                TLorentzVector bq_p4 = GetP4(genpart, b_idx);
                TLorentzVector bbarq_p4 = GetP4(genpart, bbar_idx);

                TLorentzVector bj1_p4, bj2_p4;
                if (b_jets.size() == 2)
                {
                    bj1_p4 = b_jets[0];
                    bj2_p4 = b_jets[1];
                }
                else
                {
                    continue;
                }

                TLorentzVector genmet;
                genmet.SetPtEtaPhiM(GenMET_pt, 0, GenMET_phi, 0);

                int hWW_idx = sig[SIG::H_WW];
                TLorentzVector HtoWW = GetP4(genpart, hWW_idx); 

                double def_mass = (HtoWW + bj1_p4 + bj2_p4).M();
                hm.Fill(mass_ideal_HWW, def_mass);

                auto input = { bj1_p4, bj2_p4, lq1_p4, lq2_p4, l_p4, genmet };
                // auto input = { bj1_p4, bj2_p4, lq1_p4, lq2_p4, genmet }; // lepton p4 is absent - causes bugs
                auto hme = Experimental::EstimateMassIdealHWW(input, lead_bjet_pt_pdf, rg, i, HtoWW);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.Fill(hme_ideal_hadW, mass);
                }

                hme = EstimateMass(input, lead_bjet_pt_pdf, rg, i);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.Fill(hme_mass, mass);
                }

                hme = Experimental::EstimateMassIdealNu2dPDF(input, pdf2d, rg, i, nu_p4);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.Fill(hme_mass_pdf2d_ideal_nu, mass);
                }

                hme = Experimental::EstimateMassIdealNu(input, lead_bjet_pt_pdf, rg, i, nu_p4);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.Fill(hme_mass_ideal_nu, mass);
                }
            }
        }
    }

    std::cout << "Finished processing\n";
    hm.Draw();
    hm.DrawStack({mass_ideal_HWW, hme_mass}, "Impact of b jet rescaling on X mass", "hme_vs_mass_ideal_HWW");
    hm.DrawStack({mass_ideal_HWW, hme_ideal_hadW}, "HME with ideal W->qq vs default mass", "hme_ideal_hadW_vs_mass_ideal_HWW");
    hm.DrawStack({mass_ideal_HWW, hme_mass_ideal_nu}, "HME with ideal nu vs default mass", "hme_ideal_nu_vs_mass_ideal_HWW");
    hm.DrawStack({mass_ideal_HWW, hme_mass_pdf2d_ideal_nu}, "Impact of b jet rescaling with 2D PDF on X mass", "hme_pdf2d_ideal_nu_vs_mass_ideal_HWW");
    hm.DrawStack({mass_ideal_HWW, hme_mass_pdf2d_ideal_nu, hme_mass_ideal_nu, hme_mass}, "Comparison of various methods", "cmp");

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << std::setprecision(3);
    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n";

    return 0;
}
