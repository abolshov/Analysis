#include <iostream>
#include <vector>

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

    TH1F* w_mass = new TH1F("w_mass", "W mass from lepton and corrected met", 31, 0.0, 150.0);
    // TH2F* met_vs_nu_px_corr = new TH2F("met_vs_nu_px_corr", "MET corrected px vs nu px", nbins, -200.0, 200.0, nbins, -200.0, 200.0);
    // TH2F* met_vs_nu_py_corr = new TH2F("met_vs_nu_py_corr", "MET corrected py vs nu py", nbins, -200.0, 200.0, nbins, -200.0, 200.0);
    // TH2F* met_vs_nu_px = new TH2F("met_vs_nu_px", "MET px vs nu px", nbins, -200.0, 200.0, nbins, -200.0, 200.0);
    // TH2F* met_vs_nu_py = new TH2F("met_vs_nu_py", "MET py vs nu py", nbins, -200.0, 200.0, nbins, -200.0, 200.0);

    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector b1, b2, j1, j2, l, nu, met;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);

        TLorentzVector met_corr;

        TH1F* w_from_lmet = new TH1F("w_from_lmet", "W mass from lepton and corrected met", nbins, 0.0, 1700.0);

        rg.SetSeed(0);
        for (size_t j = 0; j < nIter; ++j)
        {
            float nu_eta = rg.Uniform(-6, 6);
            // std::cout << "nu_eta = " << nu_eta << std::endl;
            met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);
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
            w_from_lmet->Fill((l+nu_corr).M());
            // std::cout << (l+nu_corr).M() << std::endl;
        }

        // met_vs_nu_px_corr->Fill(met_corr.Px(), nu.Px());
        // met_vs_nu_py_corr->Fill(met_corr.Py(), nu.Py());

        // met_vs_nu_px->Fill(met.Px(), nu.Px());
        // met_vs_nu_py->Fill(met.Py(), nu.Py());

        if (i % 1000 == 0)
        {
            std::string evt_dist_name = "iter_plots/w_mass/w_from_lmet_" + std::to_string(i);
            // save::save_1d_dist(w_from_lmet, evt_dist_name.c_str(), "W mass from lepton and corrected met");
        }
        // save::save_1d_dist(w_from_lmet, "w_from_lmet", "W mass from lepton and corrected met");
        int binmax = w_from_lmet->GetMaximumBin(); 
        float x = w_from_lmet->GetXaxis()->GetBinCenter(binmax);
        w_mass->Fill(x);
        if (i % 1000 == 0)
        {
            std::cout << "Event #" << i << ": w_mass = " << x << std::endl;
        }
        // std::cout << "w_mass = " << x << std::endl;
        delete w_from_lmet;
        // break;
    }

    // save::save_2d_dist(met_vs_nu_px_corr, "corr_met_vs_nu_px", "MET px corr [GeV]", "nu px [GeV]");
    // save::save_2d_dist(met_vs_nu_py_corr, "corr_met_vs_nu_py", "MET py corr [GeV]", "nu py [GeV]");

    // save::save_2d_dist(met_vs_nu_px, "met_vs_nu_px", "MET px [GeV]", "nu px [GeV]");
    // save::save_2d_dist(met_vs_nu_py, "met_vs_nu_py", "MET py [GeV]", "nu py [GeV]");


    save::save_1d_dist(w_mass, "w_mass", "W mass from lepton and corrected met");

    return 0;
}