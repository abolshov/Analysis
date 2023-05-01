#include <iostream>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include <math.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TSpline.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <functional>
#include "TLegend.h"
#include <numeric>
#include "TRandom3.h"
#include <cmath>
#include <random>
#include <memory>

#include "tools.hpp"


int main()
{
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
    TDirectory *dir = (TDirectory*)myFile;
    TTree *myTree = (TTree*)dir->Get("syncTree");

    float genw_mass;

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

    myTree->SetBranchAddress("genw1_mass", &genw_mass);

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

    std::unique_ptr<TH1F> onshell_w_from_qq = std::make_unique<TH1F>("onshell_w_from_qq", "On-shell W mass from gen quarks", nbins, 0.0, 100.0);
    std::unique_ptr<TH1F> offshell_w_from_qq = std::make_unique<TH1F>("offshell_w_from_qq", "Off-shell W mass from gen quarks", nbins, 0.0, 100.0);
    std::unique_ptr<TH1F> lead_on = std::make_unique<TH1F>("lead_on", "rescale PDF for leading onshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> lead_off = std::make_unique<TH1F>("lead_off", "rescale PDF for leading offshell jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> lead_bjet_pdf = std::make_unique<TH1F>("lead_bjet_pdf", "rescale PDF for leading b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> nu_eta = std::make_unique<TH1F>("nu_eta", "Neutrino rapidity distribution", nbins, -6.0, 6.0);
    std::unique_ptr<TH1F> h_mass = std::make_unique<TH1F>("h_mass", "Higgs mass", nbins, 124.92, 125.12);
    // TH1F* onshell_w_from_qq = new TH1F("onshell_w_from_qq", "On-shell W mass from gen quarks", nbins, 0.0, 100.0);
    // TH1F* offshell_w_from_qq = new TH1F("offshell_w_from_qq", "Off-shell W mass from gen quarks", nbins, 0.0, 100.0);
    // TH1F* lead_on = new TH1F("lead_on", "rescale PDF for leading onshell jet", nbins, 0.0, 6.0);
    // TH1F* lead_off = new TH1F("lead_off", "rescale PDF for leading offshell jet", nbins, 0.0, 6.0);
    // TH1F* lead_bjet_pdf = new TH1F("lead_bjet_pdf", "rescale PDF for leading b jet", nbins, 0.0, 6.0);
    // TH1F* nu_eta = new TH1F("nu_eta", "Neutrino rapidity distribution", nbins, -6.0, 6.0);
    // TH1F* h_mass = new TH1F("h_mass", "Higgs mass", nbins, 124.92, 125.12);

    int nEvents = myTree->GetEntries();
    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector q1, q2, j1, j2, l, nu, bj1, bj2, bq1, bq2;;
        q1.SetPtEtaPhiM(genq1_pt, genq1_eta, genq1_phi, genq1_mass);
        q2.SetPtEtaPhiM(genq2_pt, genq2_eta, genq2_phi, genq2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        bj1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        bj2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);

        float qq_mass = (q1 + q2).M();
        if ((l+nu).M() > qq_mass)
        {
            offshell_w_from_qq->Fill(qq_mass);
        }
        else
        {
            onshell_w_from_qq->Fill(qq_mass);
        }

        jet::pt_order(q1, q2);
        jet::pt_order(j1, j2);
        jet::pt_order(bj1, bj2);
        jet::pt_order(bq1, bq2);
        if (jet::is_offshell(q1, q2, l, nu))
        {
            lead_off->Fill(q1.Pt()/j1.Pt());
        }
        else
        {
            lead_on->Fill(q1.Pt()/j1.Pt());
        }

        lead_bjet_pdf->Fill(bq1.Pt()/bj1.Pt());
        nu_eta->Fill(gennu_eta);

        float mh = rg.Gaus(125.03, 0.03);
        h_mass->Fill(mh);
    }

    TFile *file = new TFile("pdfs.root", "RECREATE");
    file->WriteObject(lead_bjet_pdf.get(), "lead_bjet_pdf");
    file->WriteObject(offshell_w_from_qq.get(), "offshell_w_from_qq");
    file->WriteObject(onshell_w_from_qq.get(), "onshell_w_from_qq");
    file->WriteObject(lead_off.get(), "lead_off");
    file->WriteObject(lead_on.get(), "lead_on");
    file->WriteObject(h_mass.get(), "h_mass");
    file->WriteObject(nu_eta.get(), "nu_eta");
    file->Close();
    return 0;
}