#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <functional>
#include <memory>

#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TRandom3.h"

#include "Utils.hpp"
#include "Plotting.hpp"
#include "BinaryFunction.hpp"

extern TLorentzVector const zero(0.0f, 0.0f, 0.0f, 0.0f);

inline bool is_offshell(TLorentzVector const& p1, TLorentzVector const& p2, TLorentzVector const& l, TLorentzVector const& nu) noexcept
{
    return ((p1 + p2).M() < (l + nu).M());
}

void pt_order(TLorentzVector& p1, TLorentzVector& p2)
{
    if (p1.Pt() < p2.Pt())
    {
        TLorentzVector tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
}

int main()
{
    // TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-350_narrow_13TeV_Friend.root");
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_matched.root");
    // TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-800_narrow_13TeV_matched.root");
    TDirectory *dir = (TDirectory*)myFile;
    TTree *myTree = (TTree*)dir->Get("syncTree");

    float genMET_pT;
    float genMET_phi;

    float genbjet1_mass;
    float genbjet1_pt;
    float genbjet1_eta;
    float genbjet1_phi;

    float genbjet2_mass;
    float genbjet2_pt;
    float genbjet2_eta;
    float genbjet2_phi;

    float genb1_mass;
    float genb1_pt;
    float genb1_eta;
    float genb1_phi;

    float genb2_mass;
    float genb2_pt;
    float genb2_eta;
    float genb2_phi;

    float genq2_mass;
    float genq2_pt;
    float genq2_eta;
    float genq2_phi;

    float genq1_mass;
    float genq1_pt;
    float genq1_eta;
    float genq1_phi;

    float genqjet2_mass;
    float genqjet2_pt;
    float genqjet2_eta;
    float genqjet2_phi;

    float genqjet1_mass;
    float genqjet1_pt;
    float genqjet1_eta;
    float genqjet1_phi;

    float genl_mass;
    float genl_pt;
    float genl_eta;
    float genl_phi;

    float gennu_mass;
    float gennu_pt;
    float gennu_eta;
    float gennu_phi;

    myTree->SetBranchAddress("genbjet1_mass", &genbjet1_mass);
    myTree->SetBranchAddress("genbjet1_pt", &genbjet1_pt);
    myTree->SetBranchAddress("genbjet1_eta", &genbjet1_eta);
    myTree->SetBranchAddress("genbjet1_phi", &genbjet1_phi);

    myTree->SetBranchAddress("genbjet2_mass", &genbjet2_mass);
    myTree->SetBranchAddress("genbjet2_pt", &genbjet2_pt);
    myTree->SetBranchAddress("genbjet2_eta", &genbjet2_eta);
    myTree->SetBranchAddress("genbjet2_phi", &genbjet2_phi);

    myTree->SetBranchAddress("genb1_mass", &genb1_mass);
    myTree->SetBranchAddress("genb1_pt", &genb1_pt);
    myTree->SetBranchAddress("genb1_eta", &genb1_eta);
    myTree->SetBranchAddress("genb1_phi", &genb1_phi);

    myTree->SetBranchAddress("genb2_mass", &genb2_mass);
    myTree->SetBranchAddress("genb2_pt", &genb2_pt);
    myTree->SetBranchAddress("genb2_eta", &genb2_eta);
    myTree->SetBranchAddress("genb2_phi", &genb2_phi);

    myTree->SetBranchAddress("genq2_mass", &genq2_mass);
    myTree->SetBranchAddress("genq2_pt", &genq2_pt);
    myTree->SetBranchAddress("genq2_eta", &genq2_eta);
    myTree->SetBranchAddress("genq2_phi", &genq2_phi);

    myTree->SetBranchAddress("genq1_mass", &genq1_mass);
    myTree->SetBranchAddress("genq1_pt", &genq1_pt);
    myTree->SetBranchAddress("genq1_eta", &genq1_eta);
    myTree->SetBranchAddress("genq1_phi", &genq1_phi);

    myTree->SetBranchAddress("genqjet1_mass", &genqjet1_mass);
    myTree->SetBranchAddress("genqjet1_pt", &genqjet1_pt);
    myTree->SetBranchAddress("genqjet1_eta", &genqjet1_eta);
    myTree->SetBranchAddress("genqjet1_phi", &genqjet1_phi);

    myTree->SetBranchAddress("genqjet2_mass", &genqjet2_mass);
    myTree->SetBranchAddress("genqjet2_pt", &genqjet2_pt);
    myTree->SetBranchAddress("genqjet2_eta", &genqjet2_eta);
    myTree->SetBranchAddress("genqjet2_phi", &genqjet2_phi);
    
    myTree->SetBranchAddress("genl1_mass", &genl_mass);
    myTree->SetBranchAddress("genl1_pt", &genl_pt);
    myTree->SetBranchAddress("genl1_eta", &genl_eta);
    myTree->SetBranchAddress("genl1_phi", &genl_phi);

    myTree->SetBranchAddress("gennu1_mass", &gennu_mass);
    myTree->SetBranchAddress("gennu1_pt", &gennu_pt);
    myTree->SetBranchAddress("gennu1_eta", &gennu_eta);
    myTree->SetBranchAddress("gennu1_phi", &gennu_phi);

    myTree->SetBranchAddress("genMET_pT", &genMET_pT);
    myTree->SetBranchAddress("genMET_phi", &genMET_phi);

    int nbins = 101;
    TRandom3 rg;
    rg.SetSeed(0);

    int same_bquarks = 0;
    int same_bjets = 0;
    int same_light_quarks = 0;
    int same_light_jets = 0;

    // basic 1d distributions
    std::unique_ptr<TH1F> onshell_w_from_qq = std::make_unique<TH1F>("onshell_w_from_qq", "On-shell W mass from gen quarks", nbins, 0.0, 100.0);
    std::unique_ptr<TH1F> offshell_w_from_qq = std::make_unique<TH1F>("offshell_w_from_qq", "Off-shell W mass from gen quarks", nbins, 0.0, 100.0);
    std::unique_ptr<TH1F> lead_on = std::make_unique<TH1F>("lead_on", "rescale PDF for leading onshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> sublead_on = std::make_unique<TH1F>("sublead_on", "rescale PDF for subleading onshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> lead_off = std::make_unique<TH1F>("lead_off", "rescale PDF for leading offshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> sublead_off = std::make_unique<TH1F>("sublead_off", "rescale PDF for subleading offshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> lead_bjet_pdf = std::make_unique<TH1F>("lead_bjet_pdf", "rescale PDF for leading b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> sublead_bjet_pdf = std::make_unique<TH1F>("sublead_bjet_pdf", "rescale PDF for subleading b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> nu_eta = std::make_unique<TH1F>("nu_eta", "Neutrino rapidity distribution", nbins, -6.0, 6.0);
    // std::unique_ptr<TH1F> h_mass = std::make_unique<TH1F>("h_mass", "Higgs mass", nbins, 124.92, 125.12);

    std::unique_ptr<TH1F> higgs_mass = std::make_unique<TH1F>("higgs_mass", "Higgs mass from b quarks", nbins, 123.5, 126.5);
    std::unique_ptr<TH1F> qq_dR = std::make_unique<TH1F>("qq_dR", "dR between light quarks", nbins, 0.0f, 6.0f);


    std::unique_ptr<TH2F> hadrW = std::make_unique<TH2F>("hadrW", "Hadronically decaying W", nbins, 0.0f, 120.0f, nbins, 0.0f, 150.0f);
    std::unique_ptr<TH1F> jq1_dR = std::make_unique<TH1F>("jq_dR1", "dR between genjet 1 and genquark 1", nbins, 0.0, 1.0);
    std::unique_ptr<TH1F> jq2_dR = std::make_unique<TH1F>("jq_dR2", "dR between genjet 2 and genquark 2", nbins, 0.0, 1.0);

    std::unique_ptr<TH2F> jq1_eta = std::make_unique<TH2F>("jq1_eta", "Eta of Jet vs Eta of quark #1", nbins, -6.0f, 6.0f, nbins, -6.0f, 6.0f);
    std::unique_ptr<TH2F> jq2_eta = std::make_unique<TH2F>("jq2_eta", "Eta of Jet vs Eta of quark #2", nbins, -6.0f, 6.0f, nbins, -6.0f, 6.0f);
    std::unique_ptr<TH2F> jq1_phi = std::make_unique<TH2F>("jq1_phi", "Phi of Jet vs Eta of quark #1", nbins, -4.0f, 4.0f, nbins, -4.0f, 4.0f);
    std::unique_ptr<TH2F> jq2_phi = std::make_unique<TH2F>("jq2_phi", "Phi of Jet vs Eta of quark #2", nbins, -4.0f, 4.0f, nbins, -4.0f, 4.0f);
    std::unique_ptr<TH2F> jq1_E = std::make_unique<TH2F>("jq1_E", "E of Jet vs E of quark #1", nbins, 0.0f, 300.0f, nbins, 0.0f, 300.0f);
    std::unique_ptr<TH2F> jq1_Pt = std::make_unique<TH2F>("jq1_Pt", "Pt of Jet vs Pt of quark #1", nbins, 0.0f, 300.0f, nbins, 0.0f, 300.0f);
    std::unique_ptr<TH2F> jq2_E = std::make_unique<TH2F>("jq2_E", "E of Jet vs E of quark #2", nbins, 0.0f, 300.0f, nbins, 0.0f, 300.0f);
    std::unique_ptr<TH2F> jq2_Pt = std::make_unique<TH2F>("jq2_Pt", "Pt of Jet vs Pt of quark #2", nbins, 0.0f, 300.0f, nbins, 0.0f, 300.0f);

    std::unique_ptr<TH1F> jq1_dEta = std::make_unique<TH1F>("jq1_dEta", "dEta between quark and jet #1", nbins, -2.0f, 2.0f);
    std::unique_ptr<TH1F> jq1_dPhi = std::make_unique<TH1F>("jq1_dPhi", "dPhi between quark and jet #1", nbins, -1.0f, 1.0f);
    std::unique_ptr<TH1F> jq1_dE = std::make_unique<TH1F>("jq1_dE", "dE between quark and jet #1", nbins, -150.0f, 150.0f);
    std::unique_ptr<TH1F> jq1_dPt = std::make_unique<TH1F>("jq1_dPt", "dPt between quark and jet #1", nbins, -150.0f, 150.0f);
    std::unique_ptr<TH1F> jq2_dEta = std::make_unique<TH1F>("jq2_dEta", "dEta between quark and jet #2", nbins, -2.0f, 2.0f);
    std::unique_ptr<TH1F> jq2_dPhi = std::make_unique<TH1F>("jq2_dPhi", "dPhi between quark and jet #2", nbins, -1.0f, 1.0f);
    std::unique_ptr<TH1F> jq2_dE = std::make_unique<TH1F>("jq2_dE", "dE between quark and jet #2", nbins, -150.0f, 150.0f);
    std::unique_ptr<TH1F> jq2_dPt = std::make_unique<TH1F>("jq2_dPt", "dPt between quark and jet #2", nbins, -150.0f, 150.0f);

    size_t nEvents = myTree->GetEntries();
    std::cout << "nEvents = " << nEvents << "\n";
    int failedMatching = 0;
    BinaryFunction deltaR("dR", [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaR(p2); }, 0.4f);
    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector q1, q2, j1, j2, l, nu, bj1, bj2, bq1, bq2, met;
        q1.SetPtEtaPhiM(genq1_pt, genq1_eta, genq1_phi, genq1_mass);
        q2.SetPtEtaPhiM(genq2_pt, genq2_eta, genq2_phi, genq2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        bj1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        bj2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, 4.18);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, 4.18);
        met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);

        std::vector<TLorentzVector> partons = {bq1, bq2, q1, q2, l, met};
        std::vector<TLorentzVector> genjets = {bj1, bj2, j1, j2, l, met};

        // if (j1.DeltaR(q1) > 0.3f || j2.DeltaR(q2) > 0.3f) continue;

        // if (!ValidEvent(partons) || !ValidEvent(genjets)) continue;
        // if (!ValidEvent(partons)) continue;
        // if (!PartonCut(partons)) continue;
        // if (!GenJetCut(genjets)) continue;

        
        auto dEta_1 = genqjet1_eta - genq1_eta;
        auto dEta_2 = genqjet2_eta - genq2_eta;
        auto dphi_1 = abs(genqjet1_phi - genq1_phi);
        auto dPhi_1 = dphi_1 < TMath::Pi() ? dphi_1 : 2.0*TMath::Pi() - dphi_1;
        auto dphi_2 = abs(genqjet2_phi - genq2_phi);
        auto dPhi_2 = dphi_2 < TMath::Pi() ? dphi_2 : 2.0*TMath::Pi() - dphi_2;
        auto dPt_1 = genqjet1_pt - genq1_pt;
        auto dPt_2 = genqjet2_pt - genq2_pt;
        auto dE_1 = j1.E() - q1.E();
        auto dE_2 = j2.E() - q2.E();
        
        if (deltaR(j1, q1) > 0.4f || deltaR(j2, q2) > 0.4f)
        {
            // std::cout << "Event #" << i << ":\n";
            // std::cout << "q1 = ";
            // Print(q1, false); 
            // std::cout << "q2 = ";
            // Print(q2, false);
            // std::cout << "j1 = ";
            // Print(j1, false); 
            // std::cout << "j2 = ";
            // Print(j2, false); 
            // std::cout << "dR(j1, q1) = " << deltaR(j1, q1) << "\n";
            // std::cout << "dR(j2, q2) = " << deltaR(j2, q2) << "\n";
            // std::cout << "------------------------------------------------\n";
            ++failedMatching;
            continue;
        }

        if (abs(dE_1) >= 50.0f || abs(dE_2) >= 50.0f)
        {
            std::cout << "Event #" << i << ":\n";
            std::cout << "q1 = ";
            Print(q1, false); 
            std::cout << "q2 = ";
            Print(q2, false);
            std::cout << "j1 = ";
            Print(j1, false); 
            std::cout << "j2 = ";
            Print(j2, false); 
            std::cout << "dR(j1, q1) = " << deltaR(j1, q1) << "\n";
            std::cout << "dR(j2, q2) = " << deltaR(j2, q2) << "\n";
            std::cout << "dE_1 = " << dE_1 << "\n"
                      << "dPt_1 = " << dPt_1 << "\n"
                      << "dEta_1 = " << dEta_1 << "\n"
                      << "dPhi_1 = " << dPhi_1 << "\n"
                      << "dE_2 = " << dE_2 << "\n"
                      << "dPt_2 = " << dPt_2 << "\n"
                      << "dEta_2 = " << dEta_2 << "\n"
                      << "dPhi_2 = " << dPhi_2 << "\n";
            std::cout << "bq1 = ";
            Print(bq1, false);
            std::cout << "bq2 = ";
            Print(bq2, false);
            std::cout << "bj1 = ";
            Print(bj1, false);
            std::cout << "bj2 = ";
            Print(bj2, false);
            deltaR.Print({bq1, bq2, bj1, bj2, q1, q2, j1, j2});
            std::cout << "------------------------------------------------\n";
            // break;
        }

        jq1_eta->Fill(genq1_eta, genqjet1_eta);
        jq1_dEta->Fill(dEta_1);
        jq2_eta->Fill(genq2_eta, genqjet2_eta);
        jq2_dEta->Fill(dEta_2);
        jq1_phi->Fill(genq1_phi, genqjet1_phi);
        jq1_dPhi->Fill(dPhi_1);
        jq2_phi->Fill(genq2_phi, genqjet2_phi);
        jq2_dPhi->Fill(dPhi_2);
        jq1_Pt->Fill(genq1_pt, genqjet1_pt);
        jq1_dPt->Fill(dPt_1);
        jq2_Pt->Fill(genq2_pt, genqjet2_pt);
        jq2_dPt->Fill(dPt_2);
        jq1_E->Fill(q1.E(), j1.E());
        jq2_E->Fill(q2.E(), j2.E());
        jq1_dE->Fill(dE_1);
        jq2_dE->Fill(dE_2);

        qq_dR->Fill(q1.DeltaR(q2));

        // if (q1 == q2) 
        // {
        //     ++same_light_quarks;
        //     continue;
        // }
        // if (j1 == j2) 
        // {
        //     ++same_light_jets;
        //     continue;
        // }
        // if (bq1 == bq2) 
        // {
        //     ++same_bquarks;
        //     continue;    
        // }
        // if (bj1 == bj2) 
        // {
        //     ++same_bjets;
        //     continue;
        // }
        // w_from_quarks->Fill((q1 + q2).M());

        float qq_mass = (q1 + q2).M();
        float jj_mass = (j1 + j2).M();
        hadrW->Fill(qq_mass, jj_mass);
        float dr1 = j1.DeltaR(q1);
        jq1_dR->Fill(dr1);
        float dr2 = j2.DeltaR(q2);
        jq2_dR->Fill(dr2);
        float lv_mass = (l + nu).M();
        float ww_mass = (q1 + q2 + l + nu).M();
        if (lv_mass > qq_mass)
        {
            // lv is onshell
            offshell_w_from_qq->Fill(qq_mass);
        }
        if (lv_mass < qq_mass)
        {
            // jj is onshell
            onshell_w_from_qq->Fill(qq_mass);
        }

        // higgs_from_jets->Fill((bj1 + bj2).M());

        pt_order(q1, q2);
        pt_order(j1, j2);
        pt_order(bj1, bj2);
        pt_order(bq1, bq2);
        if (is_offshell(q1, q2, l, nu))
        {
            // W->qq is ofshell => fill ofshell distributions
            lead_off->Fill(q1.Pt()/j1.Pt());
            sublead_off->Fill(q2.Pt()/j2.Pt());
        }
        else
        {
            // W->qq is onshell => fill ofshell onshell
            lead_on->Fill(q1.Pt()/j1.Pt());
            sublead_on->Fill(q2.Pt()/j2.Pt());
        }

        float lead_bj_resc = bq1.Pt()/bj1.Pt();
        float subl_bj_resc = bq2.Pt()/bj2.Pt();

        higgs_mass->Fill((bq1 + bq2).M());
        // higgs_from_jets->Fill((bj1 + bj2).M());
        lead_bjet_pdf->Fill(lead_bj_resc);
        sublead_bjet_pdf->Fill(subl_bj_resc);
        nu_eta->Fill(gennu_eta);

        // float mh = rg.Gaus(125.03, 0.03);
        // h_mass->Fill(mh);
    }

    std::cout << "Failed Matching = " << failedMatching << "\n";
    std::cout << "Overlapping Jets = " << deltaR.Count() << "\n";
    std::cout << "nEvents = " << nEvents << "\n";
    // std::cout << "same_bquarks = " << same_bquarks << "\n"
    //           << "same_bjets = " << same_bjets << "\n"
    //           << "same_light_qurks = " << same_light_quarks << "\n"
    //           << "same_light_jets = " << same_light_jets << "\n";

    // TFile *file = new TFile("new_pdfs.root", "RECREATE");
    // file->WriteObject(lead_bjet_pdf.get(), "lead_bjet_pdf");
    // file->WriteObject(sublead_bjet_pdf.get(), "sublead_bjet_pdf");
    // file->WriteObject(offshell_w_from_qq.get(), "offshell_w_from_qq");
    // file->WriteObject(onshell_w_from_qq.get(), "onshell_w_from_qq");
    // file->WriteObject(lead_off.get(), "lead_off");
    // file->WriteObject(lead_on.get(), "lead_on");
    // file->WriteObject(sublead_off.get(), "sublead_off");
    // file->WriteObject(sublead_on.get(), "sublead_on");
    // // file->WriteObject(h_mass.get(), "h_mass");
    // file->WriteObject(higgs_mass.get(), "h_mass");
    // file->WriteObject(nu_eta.get(), "nu_eta");

    save_2d_dist(hadrW.get(), "hadronic_W", "diquark mass [GeV]", "dijet mass [GeV]");
    // save_1d_dist(qq_dR.get(), "qq_dR", "dR");
    save_2d_dist(jq1_eta.get(), "jq1_eta", "quark #1 eta", "jet #1 eta");
    save_2d_dist(jq2_eta.get(), "jq2_eta", "quark #2 eta", "jet #2 eta");
    save_2d_dist(jq1_phi.get(), "jq1_phi", "quark #1 phi", "jet #1 phi");
    save_2d_dist(jq2_phi.get(), "jq2_phi", "quark #2 phi", "jet #2 phi");
    save_2d_dist(jq1_E.get(), "jq1_E", "quark #1 E [GeV]", "jet #1 E [GeV]");
    save_2d_dist(jq2_E.get(), "jq2_E", "quark #2 E [GeV]", "jet #2 E [GeV]");
    save_2d_dist(jq1_Pt.get(), "jq1_Pt", "quark #1 Pt [GeV]", "jet #1 Pt [GeV]");
    save_2d_dist(jq2_Pt.get(), "jq2_Pt", "quark #2 Pt [GeV]", "jet #2 Pt [GeV]");
    save_1d_dist(jq1_dR.get(), "jq1_dR", "dR");
    save_1d_dist(jq2_dR.get(), "jq2_dR", "dR");
    save_1d_dist(jq1_dE.get(), "jq1_dE", "[GeV]");
    save_1d_dist(jq2_dE.get(), "jq2_dE", "[GeV]");
    save_1d_dist(jq1_dPhi.get(), "jq1_dPhi", "dPhi");
    save_1d_dist(jq2_dPhi.get(), "jq2_dPhi", "dPhi");
    save_1d_dist(jq1_dEta.get(), "jq1_dEta", "dEta");
    save_1d_dist(jq2_dEta.get(), "jq2_dEta", "dEta");
    save_1d_dist(jq1_dPt.get(), "jq1_dPt", "[GeV]");
    save_1d_dist(jq2_dPt.get(), "jq2_dPt", "[GeV]");
    return 0;
}