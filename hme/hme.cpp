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

#include "tools.hpp"

int main()
{
    TFile* file_pdf = new TFile("pdfs.root", "READ");
    TH1F* lead_bjet_pdf = static_cast<TH1F*>(file_pdf->Get("lead_bjet_pdf"));
    TH1F* lead_off = static_cast<TH1F*>(file_pdf->Get("lead_off"));
    TH1F* lead_on = static_cast<TH1F*>(file_pdf->Get("lead_on"));
    TH1F* offshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("offshell_w_from_qq"));
    TH1F* onshell_w_from_qq = static_cast<TH1F*>(file_pdf->Get("onshell_w_from_qq"));
    TH1F* h_mass = static_cast<TH1F*>(file_pdf->Get("h_mass"));

    TFile* file_data = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
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
    int nIter = 1000;
    float met_sigma = 25.2;
    int nbins = 101;

    TRandom3 rg;
    float mh = 125.0;

    TH1F* hh_mass_final = new TH1F("hh_mass", "HME prediction", nbins, 0.0, 1500.0);

    auto start = std::chrono::high_resolution_clock::now();

    // int seed;
    std::ofstream file;
    file.open("debug.log");

    TLorentzVector zero(0.0, 0.0, 0.0, 0.0);
    std::vector<float> repetitions;

    for (size_t i = 2; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector b1, b2, j1, j2, l, nu, met, bq1, bq2;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);

        TLorentzVector met_corr;

        TH1F* hh_mass = new TH1F(("hh_mass_" + std::to_string(i)).c_str(), ("hh_mass_" + std::to_string(i)).c_str(), 2001, 0.0, 2000.0);

        std::ofstream evt_file;
        evt_file.open("log/debug_" + std::to_string(i) + ".txt");

        if (i % 1000 == 0)
        {
            file << "Event #" << i << ":\n";

            file << "-- b jet 1 = (" << b1.Px() << ", " << b1.Py() << ", " << b1.Pz() << ", " << b1.E() << ")\n"; 
            file << "-- b jet 2 = (" << b2.Px() << ", " << b2.Py() << ", " << b2.Pz() << ", " << b2.E() << ")\n";

            file << "-- W jet 1 = (" << j1.Px() << ", " << j1.Py() << ", " << j1.Pz() << ", " << j1.E() << ")\n"; 
            file << "-- W jet 2 = (" << j2.Px() << ", " << j2.Py() << ", " << j2.Pz() << ", " << j2.E() << ")\n";

            file << "-- lepton = (" << l.Px() << ", " << l.Py() << ", " << l.Pz() << ", " << l.E() << ")\n"; 
            file << "-- nu = (" << nu.Px() << ", " << nu.Py() << ", " << nu.Pz() << ", " << nu.E() << ")\n";

            file << "-- met = (" << met.Px() << ", " << met.Py() << ", " << met.Pz() << ", " << met.E() << ")\n";

            evt_file << "-- b jet 1 = (" << b1.Px() << ", " << b1.Py() << ", " << b1.Pz() << ", " << b1.E() << ")\n"; 
            evt_file << "-- b jet 2 = (" << b2.Px() << ", " << b2.Py() << ", " << b2.Pz() << ", " << b2.E() << ")\n";

            evt_file << "-- W jet 1 = (" << j1.Px() << ", " << j1.Py() << ", " << j1.Pz() << ", " << j1.E() << ")\n"; 
            evt_file << "-- W jet 2 = (" << j2.Px() << ", " << j2.Py() << ", " << j2.Pz() << ", " << j2.E() << ")\n";

            evt_file << "-- lepton = (" << l.Px() << ", " << l.Py() << ", " << l.Pz() << ", " << l.E() << ")\n"; 
            evt_file << "-- nu = (" << nu.Px() << ", " << nu.Py() << ", " << nu.Pz() << ", " << nu.E() << ")\n";

            evt_file << "-- met = (" << met.Px() << ", " << met.Py() << ", " << met.Pz() << ", " << met.E() << ")\n";
        }
      
        rg.SetSeed(0);
        // seed = time(NULL);
        // rg.SetSeed(seed + i);
        // loop over iterations
        for (size_t j = 0; j < nIter; ++j)
        {
            // rg.SetSeed(j);
            float nu_eta = rg.Uniform(-6, 6);
            std::pair<float, float> light_jet_resc;

            if (jet::is_offshell(j1, j2, l, nu))
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_off, offshell_w_from_qq);
            }
            else
            {
                light_jet_resc = jet::compute_resc_factors(j1, j2, lead_on, onshell_w_from_qq);
            }

            std::pair<float, float> b_jet_resc = jet::compute_resc_factors(b1, b2, lead_bjet_pdf, h_mass);

            float dpx = rg.Gaus(0, met_sigma);
            float dpy = rg.Gaus(0, met_sigma);

            float c1 = b_jet_resc.first;
            float c2 = b_jet_resc.second;
            float c3 = light_jet_resc.first;
            float c4 = light_jet_resc.second;

            float met_px_corr = -(c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - (c3 - 1)*j1.Px() - (c4 - 1)*j2.Px();
            float met_py_corr = -(c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - (c3 - 1)*j1.Py() - (c4 - 1)*j2.Py();

            met_corr.SetPxPyPzE(met.Px() + dpx + met_px_corr, met.Py() + dpy + met_py_corr, met.Pz(), met.E());

            TLorentzVector nu_corr;
            nu_corr.SetPtEtaPhiM(met_corr.Pt(), nu_eta, met_corr.Phi(), 0.0);

            float tmp_hh_mass = (c1*b1 + c2*b2 + c3*j1 + c4*j2 + nu_corr + l).M();
            hh_mass->Fill(tmp_hh_mass);

            if (i % 1000 == 0)
            {
                file << "iteration #" << j << ":\n";
                file << "---- randomly generated parameter nu_eta = " << nu_eta << std::endl;

                file << "---- c1 = " << c1 << std::endl;
                file << "---- c2 = " << c2 << std::endl;
                file << "---- c3 = " << c3 << std::endl;
                file << "---- c4 = " << c4 << std::endl;

                file << "-----------------------------------------------------------------------------------------\n";

                evt_file << "iteration #" << j << ":\n";
                evt_file << "---- randomly generated parameter nu_eta = " << nu_eta << std::endl;

                evt_file << "---- c1 = " << c1 << std::endl;
                evt_file << "---- c2 = " << c2 << std::endl;
                evt_file << "---- c3 = " << c3 << std::endl;
                evt_file << "---- c4 = " << c4 << std::endl;

                evt_file << "-----------------------------------------------------------------------------------------\n";
            }
        }

        // break;
        evt_file.close();

        int binmax = hh_mass->GetMaximumBin(); 
        float evt_hh_mass = hh_mass->GetXaxis()->GetBinCenter(binmax);
        if (i % 1000 == 0)
        {
            hme_prediction.push_back(hh_mass);
            std::cout << "Event #" << i << ": hh_mass = " << evt_hh_mass << std::endl;
        }
        hh_mass_final->Fill(evt_hh_mass);
        repetitions.push_back(evt_hh_mass);
        // delete hh_mass;
    }

    int size_before = repetitions.size();
    std::vector<float>::iterator last = std::unique(repetitions.begin(), repetitions.end());
    repetitions.erase(last, repetitions.end());
    int size_after = repetitions.size();

    std::cout << "number of mass repetitions = " << size_before - size_after << std::endl;

    file.close();

    auto stop = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < hme_prediction.size(); ++i)
    {   
        std::string name = "hh_mass_" + std::to_string((i + 1)*1000);
        save::save_1d_dist(hme_prediction[i], name.c_str(), name.c_str());
    }

    auto duration = duration_cast<std::chrono::seconds>(stop - start);
 
    std::cout << "Time taken by HME event loop: " << duration.count() << " seconds" << std::endl;

    int binmax = hh_mass_final->GetMaximumBin(); 
    float hh_mass = hh_mass_final->GetXaxis()->GetBinCenter(binmax);

    std::cout << "Heavy Higgs mass = " << hh_mass << std::endl;

    save::save_1d_dist(hh_mass_final, "hh_mass", "HME prediction");
    delete file_pdf;
    return 0;
}