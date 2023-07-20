#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <time.h> 

#include "TFile.h"
#include "TDirectory.h"
#include "TTreeReader.h"
#include "TRandom3.h"
#include "TF1.h"

#include "tools.hpp"
#include "Tracer.hpp"

extern bool random_hme_failed;
extern bool simpl_impr_hme_failed;
extern float random_hme_eff;
extern float simpl_impr_hme_eff;

int main()
{
    TFile* file_pdf = new TFile("pdfs.root", "READ");
    TH1F* lead_bjet_pdf = static_cast<TH1F*>(file_pdf->Get("lead_bjet_pdf"));
    TH1F* lead_off = static_cast<TH1F*>(file_pdf->Get("lead_off"));
    TH1F* lead_on = static_cast<TH1F*>(file_pdf->Get("lead_on"));
    TH1F* offshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("offshell_w_from_qq"));
    TH1F* onshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("onshell_w_from_qq"));
    TH1F* h_mass = static_cast<TH1F*>(file_pdf->Get("h_mass"));
    TH1F* nu_eta = static_cast<TH1F*>(file_pdf->Get("nu_eta"));

    TFile* file_data = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
    // TFile* file_data = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-800_narrow_13TeV_Friend.root");
    TDirectory *dir = (TDirectory*)file_data;
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

    int nEvents = myTree->GetEntries();

    std::vector<TH1F*> hme_prediction;
    int nIter = 2000;
    float met_sigma = 25.2;
    int nbins = 101;

    TRandom3 rg;
    float mh = 125.0;

    TH1F* hh_mass_final = new TH1F("hh_mass", "Random Sampling HME", 40, 0.0, 2000.0);
    TH1F* hh_mass_final_no_light = new TH1F("hh_mass", "Random Sampling HME without light jets", 40, 0.0, 2000.0);
    TH1F* hh_mass_final_uni_nu = new TH1F("hh_mass", "Random Sampling HME with uniform nu eta", 40, 0.0, 2000.0);
    TH1F* hh_mass_final_simple = new TH1F("hh_mass", "Random Sampling HME without light jets and uniform eta", 40, 0.0, 2000.0);
    TH1F* hme_mass_final_simplified = new TH1F("hh_mass", "simplified HME", 40, 0.0, 2000.0);
    TH1F* hme_mass_final_simpl_impr = new TH1F("hh_mass", "Simplified Improved HME", 40, 0.0, 2000.0);

    auto start = std::chrono::high_resolution_clock::now();

    TLorentzVector zero(0.0, 0.0, 0.0, 0.0);
    int fails = 0;
    bool uniform_eta, light_on;

    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        std::string msg("Event #");
        msg += std::to_string(i) + ":\n";

        TLorentzVector b1, b2, j1, j2, l, nu, met, bq1, bq2;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);

        std::vector<TLorentzVector> particles = {b1, b2, j1, j2, l, nu, met};
        std::vector<TH1F*> pdfs = {lead_bjet_pdf, lead_on, lead_off, onshell_w_from_qq, offshell_w_from_qq, h_mass, nu_eta};

        rg.SetSeed(0);

        // Random Sampling HME
        light_on = true;
        uniform_eta = false;
        float evt_hh_mass = hme::hme_rand_sampl(particles, pdfs, nIter, rg, light_on, uniform_eta);
        hh_mass_final->Fill(evt_hh_mass);
        
        if (random_hme_failed)
        {
            msg += "\trandom sampling HME efficiency = " + std::to_string(random_hme_eff);
            Tracer::instance().write(msg);
        }

        // Random Sampling without light jets
        light_on = false;
        uniform_eta = false;
        float evt_hh_mass_no_light = hme::hme_rand_sampl(particles, pdfs, nIter, rg, light_on, uniform_eta);
        hh_mass_final_no_light->Fill(evt_hh_mass_no_light);
        
        if (random_hme_failed)
        {
            msg += "\trandom sampling HME without light jets efficiency = " + std::to_string(random_hme_eff);
            Tracer::instance().write(msg);
        }

        // Random Sampling with uniform eta and without light jets
        light_on = false;
        uniform_eta = true;
        float evt_hh_mass_simple = hme::hme_rand_sampl(particles, pdfs, nIter, rg, light_on, uniform_eta);
        hh_mass_final_simple->Fill(evt_hh_mass_simple);
        
        if (random_hme_failed)
        {
            msg += "\trandom sampling HME without light jets and uniform eta efficiency = " + std::to_string(random_hme_eff);
            Tracer::instance().write(msg);
        }

        // Random Sampling with uniform eta and with light jets
        light_on = true;
        uniform_eta = true;
        float evt_hh_mass_uni= hme::hme_rand_sampl(particles, pdfs, nIter, rg, light_on, uniform_eta);
        hh_mass_final_uni_nu->Fill(evt_hh_mass_uni);
        
        if (random_hme_failed)
        {
            msg += "\trandom sampling HME with light jets and uniform eta efficiency = " + std::to_string(random_hme_eff);
            Tracer::instance().write(msg);
        }

        // Simplified improved hme
        evt_hh_mass = hme::hme_simpl_impr({b1, b2, j1, j2, l, met}, h_mass, nIter, rg);
        hme_mass_final_simpl_impr->Fill(evt_hh_mass);

        // if (simpl_impr_hme_failed)
        // {
        //     msg += "\tsimplified improved HME efficiency = " + std::to_string(simpl_impr_hme_eff);
        //     Tracer::instance().write(msg);
        // }

        // simplified hme starts here
        evt_hh_mass = hme::hme_simplified({b1, b2, j1, j2, l, met});
        if (evt_hh_mass >= 0.0f) 
        {   
            // msg += "\tsimplified HME failed";
            // Tracer::instance().write(msg);
            // continue;
            hme_mass_final_simplified->Fill(evt_hh_mass);
        }
        else 
        {
            ++fails;
        }
        // hme_mass_final_simplified->Fill(evt_hh_mass);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by HME event loop: " << duration.count() << " seconds" << std::endl;

    int binmax = hh_mass_final->GetMaximumBin(); 
    float hh_mass = hh_mass_final->GetXaxis()->GetBinCenter(binmax);

    // save::save_1d_stack({hh_mass_final, hme_mass_final_simplified, hme_mass_final_simpl_impr}, 
    //                     {"random sampling HME", "simplified HME", "simplified improved HME"}, 
    //                     "hme_comparison", "HME comparison", "[GeV]");
    save::save_1d_stack({hh_mass_final, hme_mass_final_simplified}, 
                        {"random sampling HME", "simplified HME"}, 
                        "hme_comparison", "HME comparison", "[GeV]");

    // save::save_1d_stack({hh_mass_final, hme_mass_final_simplified}, 
    //                     {"random sampling HME", "simplified HME"}, 
    //                     "hme_comparison", "HME comparison", "[GeV]");

    save::save_1d_stack({hh_mass_final, hh_mass_final_no_light, hh_mass_final_uni_nu, hh_mass_final_simple}, 
                        {"lj on + non uni eta", "lj off + non uni eta", "lj on + uni eta", "lj off + uni eta"}, 
                        "rand_sampl_hme_comparison", "Random Sampling HME comparison", "[GeV]");

    std::cout << "Random Sampling HME Heavy Higgs mass = " << hh_mass << std::endl;

    auto rand_par = save::save_fit(hh_mass_final, "hh_mass_fit", "[GeV]");
    auto simpl_par = save::save_fit(hme_mass_final_simplified, "hh_mass_fit_simpl", "[GeV]");

    std::cout << "Widths ratio = " << simpl_par.second / rand_par.second << std::endl;
    // save::save_1d_dist(hh_mass_final, "hh_mass", "[GeV]");

    // binmax = hme_mass_final_simplified->GetMaximumBin(); 
    // hh_mass = hme_mass_final_simplified->GetXaxis()->GetBinCenter(binmax);
    // std::cout << "simplified HME Heavy Higgs mass = " << hh_mass << std::endl;

    // save::save_1d_dist(hme_mass_final_simplified, "hh_mass_simplified", "Simplified HME prediction");
    std::cout << "Analytical solution fails in " << fails << " out of " << nEvents << " total events" << std::endl;

    // binmax = hme_mass_final_simpl_impr->GetMaximumBin(); 
    // hh_mass = hme_mass_final_simpl_impr->GetXaxis()->GetBinCenter(binmax);
    // std::cout << "simplified improved HME Heavy Higgs mass = " << hh_mass << std::endl;
    // save::save_1d_dist(hme_mass_final_simpl_impr, "hh_mass_simpl_impr", "Simplified Improved HME prediction");

    std::cout << "done\n";

    delete file_pdf;
    return 0;
}