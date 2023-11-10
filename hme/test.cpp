#include <iostream>
#include <vector>
#include <chrono>

#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TGraph.h"

using namespace std::chrono;

float hme_simplified(std::vector<TLorentzVector> const& particles)
{
    float mass = -1.0;

    float mh = 125.0;
    TLorentzVector b1 = particles[0];
    TLorentzVector b2 = particles[1];
    TLorentzVector j1 = particles[2];
    TLorentzVector j2 = particles[3];
    TLorentzVector l = particles[4];
    TLorentzVector met = particles[5];

    TLorentzVector vis(l);
    vis += j1;
    vis += j2;
    vis += b1;
    vis += b2;
    float a = mh*mh - vis.M()*vis.M() + 2.0*vis.Px()*met.Px() + 2.0*vis.Py()*met.Py();
    float A = 4.0*(vis.E()*vis.E() - vis.Pz()*vis.Pz());
    float B = -4.0*a*vis.Pz();
    float C = -4.0*vis.E()*vis.E()*(met.Px()*met.Px() + met.Py()*met.Py()) - a*a;
    float delta = B*B - 4.0*A*C;

    if (delta < 0.0f) return mass;

    float pz_1 = (-B + sqrt(delta))/(2.0*A);
    float pz_2 = (-B - sqrt(delta))/(2.0*A);

    TLorentzVector nu1, nu2;
    nu1.SetPxPyPzE(met.Px(), met.Py(), pz_1, met.E());
    nu2.SetPxPyPzE(met.Px(), met.Py(), pz_2, met.E());

    TLorentzVector jj(j1);
    jj += j2;

    TLorentzVector h_nu1(l);
    h_nu1 += nu1;
    h_nu1 += jj;
 
    TLorentzVector h_nu2(l);
    h_nu2 += nu2;
    h_nu2 += jj;
    
    // TLorentzVector nu;
    TLorentzVector tmp_hh_mom(b1);
    tmp_hh_mom += b2;

    if (abs(h_nu1.M() - mh) < abs(h_nu2.M() - mh))
    {
        // nu = nu1;
        tmp_hh_mom += h_nu1;
    }
    else
    {
        // nu = nu2;
        tmp_hh_mom += h_nu2;
    }

    double tmp_hh_mass = tmp_hh_mom.M();
    mass = tmp_hh_mass;
    return mass;
}

int main()
{
    auto prog_start = std::chrono::high_resolution_clock::now();
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J_Friend.root");
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

    size_t nEvents = myTree->GetEntries();
    auto loop_start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector b1, b2, j1, j2, l, nu, met, q1, q2;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        met.SetPtEtaPhiM(genMET_pT, 0.0, genMET_phi, 0.0);

        hme_simplified({b1, b2, j1, j2, l, met});
    }
    auto loop_stop = std::chrono::high_resolution_clock::now();
    auto loop_duration = duration_cast<std::chrono::milliseconds>(loop_stop - loop_start);
    auto prog_duration = duration_cast<std::chrono::milliseconds>(loop_stop - prog_start);
    std::cout << "Time taken by simple HME event loop: " << loop_duration.count() << " ms\n";
    std::cout << "Time taken by entire program: " << prog_duration.count() << " ms\n";


    return 0;
}