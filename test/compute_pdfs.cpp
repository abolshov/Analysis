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
#include <fstream>

#include "tools.hpp"
#include "mylib.hpp"


int main()
{
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-800_narrow_13TeV_matched.root");
    TDirectory *dir = (TDirectory*)myFile;
    TTree *myTree = (TTree*)dir->Get("syncTree");

    // int genb1_statusFlags;
    // int genb2_statusFlags;

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

    // myTree->SetBranchAddress("genb1_statusFlags", &genb1_statusFlags);
    // myTree->SetBranchAddress("genb2_statusFlags", &genb2_statusFlags);

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

    std::unique_ptr<TH1F> bqbj_lead_dR = std::make_unique<TH1F>("bqbj_lead_dR", "dR between leading b quark and b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> bqbj_subl_dR = std::make_unique<TH1F>("bqbj_subl_dR", "dR between subleading b quark and b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> bad_dR_1 = std::make_unique<TH1F>("bad_dR_1", "dR between b lead quark and b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> bad_dR_2 = std::make_unique<TH1F>("bad_dR_2", "dR between b subl quark and b jet", nbins, 0.0, 6.0);
    std::unique_ptr<TH1F> very_bad_dR_1 = std::make_unique<TH1F>("very_bad_dR_1", "dR between b lead quark and b jet", 31, 0.4, 6.0);
    std::unique_ptr<TH1F> very_bad_dR_2 = std::make_unique<TH1F>("very_bad_dR_2", "dR between b subl quark and b jet", 31, 0.4, 6.0);
    std::unique_ptr<TH1F> higgs_mass = std::make_unique<TH1F>("higgs_mass", "Higgs mass from b quarks", nbins, 123.5, 126.5);

    int nEvents = myTree->GetEntries();
    // int bad_lead_bj_resc = 0;
    // int above_dR_thresh_lead = 0;
    // int bad_subl_bj_resc = 0;
    // int above_dR_thresh_subl = 0;
    // int failed_matching = 0;

    // int bq1_lastCopyBeforeFSR, bq1_lastCopy, bq1_firstCopy;
    // int bq2_lastCopyBeforeFSR, bq2_lastCopy, bq2_firstCopy;
    // bq1_lastCopyBeforeFSR = bq1_lastCopy = bq1_firstCopy = 0;
    // bq2_lastCopyBeforeFSR = bq2_lastCopy = bq2_firstCopy = 0;

    // std::ofstream subl_file;
    // subl_file.open("subl_log.txt");
    // std::ofstream lead_file;
    // lead_file.open("lead_log.txt");

    bool low_pt, med_pt, high_pt;

    size_t count = 0;
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
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, 4.18);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, 4.18);

        float qq_mass = (q1 + q2).M();
        float lv_mass = (l + nu).M();
        float ww_mass = (q1 + q2 + l + nu).M();
        if (lv_mass > qq_mass)
        {
            offshell_w_from_qq->Fill(qq_mass);
        }
        if (lv_mass < qq_mass)
        {
            onshell_w_from_qq->Fill(qq_mass);
        }

        //onshell_w_from_qq->Fit("gaus");

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

        float lead_bj_resc = bq1.Pt()/bj1.Pt();
        float subl_bj_resc = bq2.Pt()/bj2.Pt();
        // float lead_dR = bq1.DeltaR(bj1);
        // float subl_dR = bq2.DeltaR(bj2);
        // if (lead_dR > 0.4) bqbj_lead_dR->Fill(lead_dR);
        // if (subl_dR > 0.4) bqbj_subl_dR->Fill(subl_dR);

        // if (bits::extract_bit(genb2_statusFlags, 14)) ++bq2_lastCopyBeforeFSR;
        // if (bits::extract_bit(genb2_statusFlags, 13)) ++bq2_lastCopy;
        // if (bits::extract_bit(genb2_statusFlags, 12)) ++bq2_firstCopy;

        // if (bits::extract_bit(genb1_statusFlags, 14)) ++bq1_lastCopyBeforeFSR;
        // if (bits::extract_bit(genb1_statusFlags, 13)) ++bq1_lastCopy;
        // if (bits::extract_bit(genb1_statusFlags, 12)) ++bq1_firstCopy;

        // if (lead_dR > 0.4 && subl_dR > 0.4) ++failed_matching; 

        // if (lead_dR < 0.4 && subl_dR < 0.4) higgs_mass->Fill((bq1 + bq2).M());
        higgs_mass->Fill((bq1 + bq2).M());
        // if (bits::extract_bit(genb2_statusFlags, 14) && bits::extract_bit(genb1_statusFlags, 14)) higgs_mass->Fill((bq1 + bq2).M());
        // if (lead_bj_resc >= 2.0f)
        // {
        //     lead_file << "Event #" << i << ": leading resc factor exceeds 2.0" << std::endl;
        //     lead_file << "b quark = (" << bq1.Px() << ", " << bq1.Py() << ", " << bq1.Pz() << "," << bq1.E() << ")" << std::endl;
        //     lead_file << "b jet = (" << bj1.Px() << ", " << bj1.Py() << ", " << bj1.Pz() << "," << bj1.E() << ")" << std::endl;
        //     float dR = bq1.DeltaR(bj1);
        //     if (dR >= 0.4f)
        //     {
        //         ++above_dR_thresh_lead;
        //     }
        //     lead_file << "dR = " << dR << std::endl;
        //     lead_file << "-----------------------------------------------\n";
        //     ++bad_lead_bj_resc;
        // }
        // if (subl_bj_resc >= 2.0f)
        // {
        //     subl_file << "Event #" << i << ": subl resc factor exceeds 2.0" << std::endl;
        //     subl_file << "bq1 = (" << bq1.Px() << ", " << bq1.Py() << ", " << bq1.Pz() << ", " << bq1.E() << ")" << std::endl;
        //     subl_file << "bj1 = (" << bj1.Px() << ", " << bj1.Py() << ", " << bj1.Pz() << ", " << bj1.E() << ")" << std::endl;
        //     subl_file << "bq2 = (" << bq2.Px() << ", " << bq2.Py() << ", " << bq2.Pz() << ", " << bq2.E() << ")" << std::endl;
        //     subl_file << "bj2 = (" << bj2.Px() << ", " << bj2.Py() << ", " << bj2.Pz() << ", " << bj2.E() << ")" << std::endl;
        //     float dR_2 = bq2.DeltaR(bj2);
        //     float dR_1 = bq1.DeltaR(bj1);
        //     bad_dR_1->Fill(dR_1);
        //     bad_dR_2->Fill(dR_2);
        //     if (/*dR_1 >= 0.4f || */dR_2 >= 0.4f)
        //     {
        //         ++above_dR_thresh_subl;
        //         very_bad_dR_1->Fill(dR_1);
        //         very_bad_dR_2->Fill(dR_2);
        //     }
        //     subl_file << "dR_1 = " << dR_1 << std::endl;
        //     subl_file << "dR_2 = " << dR_2 << std::endl;

        //     subl_file << "bq2: isLastCopyBeforeFSR = " << bits::extract_bit(genb2_statusFlags, 14) << std::endl;
        //     // if (bits::extract_bit(genb2_statusFlags, 14)) ++bq2_lastCopyBeforeFSR;
        //     subl_file << "bq2: isLastCopy = " << bits::extract_bit(genb2_statusFlags, 13) << std::endl;
        //     // if (bits::extract_bit(genb2_statusFlags, 13)) ++bq2_lastCopy;
        //     subl_file << "bq2: isFirstCopy = " << bits::extract_bit(genb2_statusFlags, 12) << std::endl;
        //     // if (bits::extract_bit(genb2_statusFlags, 12)) ++bq2_firstCopy;

        //     subl_file << "bq1: isLastCopyBeforeFSR = " << bits::extract_bit(genb1_statusFlags, 14) << std::endl;
        //     // if (bits::extract_bit(genb1_statusFlags, 14)) ++bq1_lastCopyBeforeFSR;
        //     subl_file << "bq1: isLastCopy = " << bits::extract_bit(genb1_statusFlags, 13) << std::endl;
        //     // if (bits::extract_bit(genb1_statusFlags, 13)) ++bq1_lastCopy;
        //     subl_file << "bq1: isFirstCopy = " << bits::extract_bit(genb1_statusFlags, 12) << std::endl;
        //     // if (bits::extract_bit(genb1_statusFlags, 12)) ++bq1_firstCopy;

        //     subl_file << "-----------------------------------------------\n";
        //     ++bad_subl_bj_resc;
        // }
        lead_bjet_pdf->Fill(lead_bj_resc);
        nu_eta->Fill(gennu_eta);

        float mh = rg.Gaus(125.03, 0.03);
        h_mass->Fill(mh);
    }

    // std::cout << "bad_lead_bj_resc = " << bad_lead_bj_resc << std::endl;
    // std::cout << "above_dR_thresh_lead = " << above_dR_thresh_lead << std::endl;
    // std::cout << "bad_subl_bj_resc = " << bad_subl_bj_resc << std::endl;
    // std::cout << "above_dR_thresh_subl = " << above_dR_thresh_subl << std::endl;
    // std::cout << "failed_matching = " << failed_matching << std::endl;
    // std::cout << "--------------------------------\n";
    // std::cout << "bq1_firstCopy = " << bq1_firstCopy << std::endl;
    // std::cout << "bq1_lastCopy = " << bq1_lastCopy << std::endl;
    // std::cout << "bq1_lastCopyBeforeFSR = " << bq1_lastCopyBeforeFSR << std::endl;
    // std::cout << "bq2_firstCopy = " << bq2_firstCopy << std::endl;
    // std::cout << "bq2_lastCopy = " << bq2_lastCopy << std::endl;
    // std::cout << "bq2_lastCopyBeforeFSR = " << bq2_lastCopyBeforeFSR << std::endl;

    save::save_1d_dist(nu_eta.get(), "nu_eta", "Eta");
    save::save_1d_dist(higgs_mass.get(), "higgs_mass", "[GeV]");

    // save::save_1d_dist(bqbj_lead_dR.get(), "bqbj_lead_dR", "dR");
    // save::save_1d_dist(bqbj_subl_dR.get(), "bqbj_subl_dR", "dR");

    // save::save_1d_dist(bad_dR_1.get(), "bad_dR_1", "dR");
    // save::save_1d_dist(bad_dR_2.get(), "bad_dR_2", "dR");

    // save::save_1d_dist(very_bad_dR_1.get(), "very_bad_dR_1", "dR");
    // save::save_1d_dist(very_bad_dR_2.get(), "very_bad_dR_2", "dR");

    // save::save_1d_stack({bad_dR_1.get(), bad_dR_2.get()}, {"q1j1", "q2j2"}, "dR_comparison", "dR_comparison", "dR");
    // save::save_1d_stack({bqbj_lead_dR.get(), bqbj_subl_dR.get()}, {"leading pair dR", "subleading pair dR"}, "dR_comparison_full", "dR_comparison_full", "dR");
    // save::save_1d_stack({very_bad_dR_1.get(), very_bad_dR_2.get()}, {"q1j1", "q2j2"}, "dR_comparison_vb", "dR_comparison_vb", "dR");

    auto g1 = new TF1("m1", "crystalball", 0.0, 100);
    g1->SetParameter(1, 80);
    g1->SetParameter(2, 5);
    g1->SetParameter(3, 1);
    g1->SetParameter(4, 5);
    g1->SetParameter(5, 5);
    g1->SetParameter(0, onshell_w_from_qq->GetEntries());

    std::cout << "blah\n";

    onshell_w_from_qq->Fit(g1, "R");

    // auto params = save::save_fit(onshell_w_from_qq.get(), "onshell_w_from_qq", "[GeV]");

    // std::cout << "# Anomalies = " << count << std::endl;

    TFile *file = new TFile("pdfs_2016.root", "RECREATE");
    file->WriteObject(lead_bjet_pdf.get(), "lead_bjet_pdf");
    file->WriteObject(offshell_w_from_qq.get(), "offshell_w_from_qq");
    file->WriteObject(onshell_w_from_qq.get(), "onshell_w_from_qq");
    file->WriteObject(lead_off.get(), "lead_off");
    file->WriteObject(lead_on.get(), "lead_on");
    // file->WriteObject(h_mass.get(), "h_mass");
    file->WriteObject(higgs_mass.get(), "h_mass");
    file->WriteObject(nu_eta.get(), "nu_eta");
    file->Close();

    // subl_file.close();
    // lead_file.close();

    return 0;
}