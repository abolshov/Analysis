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
#include "constants.hpp"

using namespace std::chrono;

int main()
{
    // TFile* file_pdf = new TFile("data/pdfs.root", "READ");
    TFile* file_pdf = new TFile("data/new_pdfs.root", "READ");
    TH1F* lead_bjet_pdf = static_cast<TH1F*>(file_pdf->Get("lead_bjet_pdf"));
    TH1F* lead_off = static_cast<TH1F*>(file_pdf->Get("lead_off"));
    TH1F* lead_on = static_cast<TH1F*>(file_pdf->Get("lead_on"));
    TH1F* offshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("offshell_w_from_qq"));
    TH1F* onshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("onshell_w_from_qq"));
    TH1F* h_mass = static_cast<TH1F*>(file_pdf->Get("h_mass"));
    TH1F* nu_eta = static_cast<TH1F*>(file_pdf->Get("nu_eta"));
    // TH1F* sublead_bjet_pdf = static_cast<TH1F*>(file_pdf->Get("sublead_bjet_pdf"));
    TH1F* sublead_on = static_cast<TH1F*>(file_pdf->Get("sublead_on"));
    TH1F* sublead_off = static_cast<TH1F*>(file_pdf->Get("sublead_off"));

    // TFile* file_data = TFile::Open("data/NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
    // TFile* file_data = TFile::Open("data/GluGluToRadionToHHTo2B2WToLNu2J_M-350_narrow_13TeV_Friend.root");
    TFile* file_data = TFile::Open("data/NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_matched.root");
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

    int nEvents = static_cast<int>(myTree->GetEntries());

    int nIter = 10000;
    // int nbins = 101;
    int zero_part_event = 0;
    int identical_pair_event = 0;
    int too_low_hh_mass = 0;
    int hme_error = 0;
    int bad_light_jet_match = 0;
    int bad_b_jet_match = 0;
    int valid_event = 0;

    TRandom3 rg;

    TH1F* hh_mass_rand_sampl = new TH1F("hh_mass", "Random Sampling", 200, 0.0, 2000.0);
    TH1F* hh_mass_rand_sampl_lj = new TH1F("hh_mass", "Random Sampling with light jets", 200, 0.0, 2000.0);
    TH1F* hh_mass_analytical = new TH1F("hh_mass", "Analytical solution", 200, 0.0, 2000.0);

    TH1F* hadronic_w = new TH1F("hadronic_w", "Hadronic W mass", 51, 0.0, 120.0);
    TH1F* dR1 = new TH1F("dR1", "dR1", 51, 0.0, 1.0);
    TH1F* dR2 = new TH1F("dR2", "dR2", 51, 0.0, 1.0);

    auto start = std::chrono::high_resolution_clock::now();
    // TLorentzVector zero(0.0, 0.0, 0.0, 0.0);
    int analytical_fails = 0;
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        // std::string msg("Event #");
        // msg += std::to_string(i) + ":\n";

        TLorentzVector bq1, bq2, b1, b2, j1, j2, l, nu, met, q1, q2;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);
        q1.SetPtEtaPhiM(genq1_pt, genq1_eta, genq1_phi, genq1_mass);
        q2.SetPtEtaPhiM(genq2_pt, genq2_eta, genq2_phi, genq2_mass);
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);

        std::vector<TLorentzVector> particles = {b1, b2, j1, j2, l, nu, met};
        // std::vector<TLorentzVector> particles = {b1, b2, q1, q2, l, nu, met}; // quarks instead of jets
        std::vector<TH1F*> pdfs = {lead_bjet_pdf, lead_on, lead_off, onshell_w_from_qq, offshell_w_from_qq, h_mass, nu_eta};

        if (!HasZeroParticle(particles))
        {
            ++zero_part_event;
            // std::cout << "Skipping event #" << i << ": contains zero particle(s)" << "\n";
            continue;
        }

        if (HasIdenticalPair(j1, j2) || HasIdenticalPair(b1, b2))
        {
            ++identical_pair_event;
            // std::cout << "Skipping event #" << i << ": contains identical pair(s)" << "\n";
            continue;
        }

        if (!ValidDeltaR(q1, j1) || !ValidDeltaR(q2, j2))
        {
            ++bad_light_jet_match;
            continue;
        }

        if (!ValidDeltaR(bq1, b1) || !ValidDeltaR(bq2, b2))
        {
            ++bad_b_jet_match;
            continue;
        }

        ++valid_event;

        hadronic_w->Fill((j1 + j2).M());
        dR1->Fill(j1.DeltaR(q1));
        dR2->Fill(j2.DeltaR(q2));

        rg.SetSeed(0);
        // Ideal HME (when I know for sure what W I am using) without light jet corrections
        float evt_hh_mass = hme::rand_sampl(particles, pdfs, nIter, rg, 3000);
        if (evt_hh_mass > 0.0f)
        {
            if (evt_hh_mass < 250.0f) 
            {
                ++too_low_hh_mass;
                continue;
            }
            hh_mass_rand_sampl->Fill(evt_hh_mass);
            // std::cout << "Event " << i << ", hh_mass_rand_sampl = " << evt_hh_mass << "\n";
        }
        else 
        {
            ++hme_error;
        }

        // random sampling with usual number of bins
        evt_hh_mass = hme::rand_sampl(particles, pdfs, nIter, rg, 3000, false, true);
        if (evt_hh_mass > 0.0f)
        {
            if (evt_hh_mass < 250.0f) continue;
            hh_mass_rand_sampl_lj->Fill(evt_hh_mass);
        }
        
        // simplified hme
        evt_hh_mass = hme::analytical({b1, b2, j1, j2, l, met});
        if (evt_hh_mass > 0.0f) 
        {   
            hh_mass_analytical->Fill(evt_hh_mass);
        }
        else 
        {
            ++analytical_fails;
        }
        // break;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by HME event loop: " << duration.count() << " seconds" << "\n";

    std::cout << "Total " << nEvents << " events, from them:\n";
    std::cout << "\t have zero particle(s) " << zero_part_event << "(" << (1.0*zero_part_event)/nEvents*100 << "%)\n";
    std::cout << "\t have identical pair(s) " << identical_pair_event << "(" << (1.0*identical_pair_event)/nEvents*100 << "%)\n";
    std::cout << "\t have bad light jet match " << bad_light_jet_match << "(" << (1.0*bad_light_jet_match)/nEvents*100 << "%)\n";
    std::cout << "\t have bad b jet match " << bad_b_jet_match << "(" << (1.0*bad_b_jet_match)/nEvents*100 << "%)\n";

    std::cout << "\n";

    std::cout << "Total " << valid_event << " valid events, from them:\n";
    std::cout << "\t have hh_mass < 2*mh " << too_low_hh_mass << "(" << (1.0*too_low_hh_mass)/valid_event*100 << "%)\n";
    std::cout << "\t have HME error (hh_mass < 0) " << hme_error << "(" << (1.0*hme_error)/valid_event*100 << "%)\n";
    std::cout << "\t have analytical solution fail " << analytical_fails << "(" << (1.0*analytical_fails)/valid_event*100 << "%)\n";

    std::cout << "Number of HME predictions: " << hh_mass_rand_sampl->GetEntries() << "\n";
    std::cout << "Number of analytical solutions: " << hh_mass_analytical->GetEntries() << "\n";

    save::save_1d_stack({hh_mass_analytical, hh_mass_rand_sampl}, 
                        {"analytical solution", "random sampling"}, 
                        "rs_vs_anal", "Random sampling vs Analytical solution", "[GeV]");

    save::save_1d_stack({hh_mass_analytical, hh_mass_rand_sampl, hh_mass_rand_sampl_lj}, 
                        {"analytical solution", "random sampling", "weighted random sampling"}, 
                        "hme_comparison", "HME comparison", "[GeV]");

    auto rand_sampl_par = save::save_fit(hh_mass_rand_sampl, "hh_mass_rand_sampl_fit", "[GeV]", "landau");
    auto simpl_par = save::save_fit(hh_mass_analytical, "hh_mass_analytical_fit", "[GeV]", "landau");
    std::cout << "Widths ratio = " << simpl_par.second / rand_sampl_par.second << "\n";

    save::save_1d_dist(hh_mass_rand_sampl, "hh_mass_rand_sampl", "[GeV]");
    save::save_1d_dist(hh_mass_analytical, "hh_mass_analytical", "[GeV]");

    save::save_1d_dist(hadronic_w, "hadronic_w", "[GeV]");
    save::save_1d_dist(dR1, "dR1", "[dR]");
    save::save_1d_dist(dR2, "dR2", "[dR]");

    // std::cout << "Analytical solution fails in " << analytical_fails << " out of " << nEvents << " total events (" << (1.0*analytical_fails)/nEvents*100 << "%)" << "\n";
    std::cout << "done\n";

    delete file_pdf;
    return 0;
}