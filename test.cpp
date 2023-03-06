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
#include "TFile.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TStyle.h"
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

float const W_threshold = 55.0;
float const off_shell_mass = 40.0;
float const on_shell_mass = 80.433;
float const w_width = 2.085;

static int fails = 0;
static int total = 0;

using TwoObjFunc = std::function<float(TLorentzVector const&, TLorentzVector const&)>;

bool is_offshell(TLorentzVector const& p1, TLorentzVector const& p2, float mass)
{   
    return (p1 + p2).M() < mass;
}

bool is_offshell(TLorentzVector const& p1,
                 TLorentzVector const& p2,
                 TLorentzVector const& l,
                 TLorentzVector const& nu)
{
    return ((p1 + p2).M() < (l + nu).M());
}

// assigns object with greater pt to p1 and the other one to p2
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

// this function always assigns leading jet to p1, trailing - p2
// in the returned pair first value is correction for leading jet p1, second - for trailing jet p2
std::pair<float, float> compute_jet_resc_factors(TLorentzVector& p1, 
                                                 TLorentzVector& p2,
                                                 TH1F* lead_jet_pdf, 
                                                 TH1F* mass_pdf)
{
    pt_order(p1, p2);

    float mass = mass_pdf->GetRandom(); 
    // std::cout << "generated mass = " << mass << std::endl;
    
    float lead_rescale_factor = lead_jet_pdf->GetRandom();
    float x1 = p2.M2();
    float x2 = 2*lead_rescale_factor*(p1*p2);
    float x3 = lead_rescale_factor*lead_rescale_factor*p1.M2() - mass*mass;

    float universal_corr = mass/(p1+p2).M();

    if (x2 < 0)
    {
        std::cout << "No solution: negative scalar product of jet momenta" << std::endl; 
        return std::make_pair(universal_corr, universal_corr);
    }
    if ((x2*x2 - 4*x1*x3) < 0 || x1 == 0)
    {
        std::cout << "No solution: negative discriminiant" << std::endl;
        return std::make_pair(universal_corr, universal_corr);
    }
    float trail_rescale_factor = (-x2 + sqrt(x2*x2 - 4*x1*x3))/(2*x1);
    if (trail_rescale_factor < 0)
    {
        std::cout << "No solution: negative trailing jet rescaling factor" << std::endl;
        return std::make_pair(universal_corr, universal_corr);
    }
    return std::make_pair(lead_rescale_factor, trail_rescale_factor);
}

std::pair<float, float> compute_jet_resc_factors(TLorentzVector& p1, 
                                                 TLorentzVector& p2,
                                                 TH1F* lead_jet_pdf, 
                                                 float mass)
{
    pt_order(p1, p2);
    // std::cout << "generated mass = " << mass << std::endl;
    
    float lead_rescale_factor = lead_jet_pdf->GetRandom();
    float x1 = p2.M2();
    float x2 = 2*lead_rescale_factor*(p1*p2);
    float x3 = lead_rescale_factor*lead_rescale_factor*p1.M2() - mass*mass;

    float universal_corr = mass/(p1+p2).M();

    if (x2 < 0)
    {
        std::cout << "No solution: negative scalar product of jet momenta" << std::endl; 
        return std::make_pair(universal_corr, universal_corr);
    }
    if ((x2*x2 - 4*x1*x3) < 0 || x1 == 0)
    {
        std::cout << "No solution: negative discriminiant" << std::endl;
        return std::make_pair(universal_corr, universal_corr);
    }
    float trail_rescale_factor = (-x2 + sqrt(x2*x2 - 4*x1*x3))/(2*x1);
    if (trail_rescale_factor < 0)
    {
        std::cout << "No solution: negative trailing jet rescaling factor" << std::endl;
        return std::make_pair(universal_corr, universal_corr);
    }
    return std::make_pair(lead_rescale_factor, trail_rescale_factor);
}


TLorentzVector compute_met_corr(TLorentzVector& met,
                                TLorentzVector const& b1,
                                TLorentzVector const& b2,
                                TLorentzVector const& j1,
                                TLorentzVector const& j2,
                                std::pair<float, float> b_jet_resc_factors, 
                                std::pair<float, float> w_jet_resc_factors, 
                                TRandom3& rg, 
                                float met_sigma)
{
    float c1 = b_jet_resc_factors.first;
    float c2 = b_jet_resc_factors.second;
    float c3 = w_jet_resc_factors.first;
    float c4 = w_jet_resc_factors.second;

    float dpx = rg.Gaus(0, met_sigma);
    float dpy = rg.Gaus(0, met_sigma);

    float met_px_corr = -(c1 - 1)*b1.Px() - (c2 - 1)*b2.Px() - (c3 - 1)*j1.Px() - (c2 - 1)*j2.Px();
    float met_py_corr = -(c1 - 1)*b1.Py() - (c2 - 1)*b2.Py() - (c3 - 1)*j1.Py() - (c2 - 1)*j2.Py();

    TLorentzVector result;
    result.SetPxPyPzE(met.Px() + dpx + met_px_corr, met.Py() + dpy + met_py_corr, met.Pz(), met.E());
    return result;
}

void rescale_jets(TLorentzVector& p1, 
                  TLorentzVector& p2,
                  TH1F* lead_jet_pdf, 
                  TH1F* mass_pdf)
{
    // std::cout << "mass before rescaling = " << (p1+p2).M() << std::endl;
    auto resc_factors = compute_jet_resc_factors(p1, p2, lead_jet_pdf, mass_pdf);
    if (resc_factors.first == -1)
    {
        ++fails;
        return;
    }
    ++total; 
    p1 *= resc_factors.first;
    p2 *= resc_factors.second;
    // std::cout << "c1 = " << resc_factors.first << ", c2 = " << resc_factors.second << std::endl;
    // std::cout << "mass after rescaling = " << (p1+p2).M() << std::endl;
    // std::cout << "-------------------------------------" << std::endl;
}

TLorentzVector generate_nu(TLorentzVector const& met /*, TRandom3& rg */)
{
    TLorentzVector nu;
    // float eta = rg.Uniform(-6, 6);
    float eta = met.Eta();
    float pt = met.Pt();
    float phi = met.Phi();
    nu.SetPtEtaPhiM(pt, phi, eta, 0.0);
    return nu;
}

void fill_2obj_dist(TLorentzVector const& p1, 
                    TLorentzVector const& p2,
                    TwoObjFunc& func,
                    TH1F* dist)
{   
    float val = func(p1, p2);
    dist->Fill(val);
}

void save_1d_dist(TH1F* dist, 
                  std::string const& name,
                  std::string const& title)
{
    TCanvas* c1 = new TCanvas("c1", "c1");

    dist->GetXaxis()->SetTitle(title.c_str());
    dist->SetLineWidth(3);
    dist->DrawNormalized();
    c1->SaveAs((name + ".png").c_str());

    delete c1;
}

void save_2d_dist(TH2F* dist, 
                  std::string const& name,
                  std::string const& title_x,
                  std::string const& title_y)
{
    TCanvas* c1 = new TCanvas("c1", "c1");
    TStyle* gStyle = new TStyle();
    gStyle->SetPalette(kRainBow);

    dist->GetXaxis()->SetTitle(title_x.c_str());
    dist->GetYaxis()->SetTitle(title_y.c_str());
    dist->SetLineWidth(3);
    dist->Draw("colz");
    c1->SaveAs((name + ".png").c_str());

    delete gStyle;
    delete c1;
}

void save_1d_stack(std::vector<TH1F*> const& distrs,
                   std::vector<std::string> const& legends,
                   std::string const& name,
                   std::string const& title,
                   std::string const& axis_label)
{
    if (distrs.size() != legends.size())
    {
        std::cout << "number of legends and histograms do not match!" << std::endl;
        return;
    }
    TCanvas* c1 = new TCanvas("c1", "c1");

    THStack* stack = new THStack("stack", title.c_str());
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (size_t i = 0; i < distrs.size(); ++i)
    {
        distrs[i]->SetLineWidth(3);
        distrs[i]->SetLineColor(i + 1);
        stack->Add(distrs[i]);
        legend->AddEntry(distrs[i], legends[i].c_str(), "l");
    }

    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle(axis_label.c_str());
    legend->Draw();
    c1->SaveAs((name + ".png").c_str());

    delete stack;
    delete legend;
    delete c1;
}

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
    
    TH1F* lead_bjet_pdf = new TH1F("lead_pdf", "rescale PDF for leading b jet", nbins, 0.0, 6.0);
    TH1F* h_from_bb_corr = new TH1F("h_from_bb_corr", "H mass from corrected b jets", nbins, 124.8, 125.2);
    TH1F* h_from_bb = new TH1F("h_from_bb", "H mass from b jets", nbins, 0.0, 150.0);
    TH1F* h_from_bq = new TH1F("h_from_bq", "H mass from b quarks", nbins, 122.0, 126.0);
    TH1F* h_from_qqlv = new TH1F("h_from_qqlv", "H mass from qqlv", nbins, 0.0, 350.0);

    // TH1F* qj1qj2_dR = new TH1F("qj1qj2_dR", "dR between light jets", nbins, -5.0, 5.0);
    TH1F* bjbq_dR = new TH1F("bjbq_dR", "dR between leading b jet and b quark", nbins, -1.0, 5.0);
    TH1F* lead_on = new TH1F("lead_on", "rescale PDF for leading onshell jet", nbins, 0.0, 6.0);
    TH1F* lead_off = new TH1F("lead_off", "rescale PDF for leading offshell jet", nbins, 0.0, 6.0);

    TH1F* w_from_qq = new TH1F("w_from_qq", "W mass from gen quarks", nbins, 0.0, 120.0);
    TH1F* w_from_jj = new TH1F("w_from_jj", "W mass from gen jets", nbins, 0.0, 120.0);
    TH1F* w_from_jj_corr = new TH1F("w_from_jj_corr", "W mass from corrected gen jets", nbins, 0.0, 120.0);
    TH1F* w_from_lv = new TH1F("w_from_lv", "W mass from gen lepton and neutrino", nbins, 0.0, 120.0);

    TH1F* lead_j_pt = new TH1F("lead_j_pt", "Leading jet p_t", nbins, 0.0, 350.0);
    TH1F* lead_j_pt_corr = new TH1F("lead_j_pt_corr", "Leading b jet p_t corrected", nbins, 0.0, 350.0);
    TH1F* lead_q_pt = new TH1F("lead_q_pt", "Leading b quark p_t", nbins, 0.0, 350.0);
    THStack* pt = new THStack("pt", "b object Pt before and after correcting");
    auto pt_legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    TH1F* onshell_w_from_qq = new TH1F("onshell_w_from_qq", "On-shell W mass from gen quarks", nbins, 0.0, 100.0);
    TH1F* offshell_w_from_qq = new TH1F("offshell_w_from_qq", "Off-shell W mass from gen quarks", nbins, 0.0, 100.0);
    TH1F* onshell_w_from_jj = new TH1F("onshell_w_from_jj", "On-shell W mass from gen jets", nbins, 0.0, 120.0);
    TH1F* offshell_w_from_jj = new TH1F("offshell_w_from_jj", "Off-shell W mass from gen jets", nbins, 0.0, 120.0);

    TH2F* genW_vs_lv = new TH2F("genW_vs_lv", "Gen W mass vs W from lv", nbins, 0.0, 120.0, nbins, 0.0, 120.0);
    TH2F* genW_vs_qq = new TH2F("genW_vs_qq", "Gen W mass vs W from qq", nbins, 0.0, 120.0, nbins, 0.0, 120.0);
    TH2F* heavy_vs_light = new TH2F("heavy_vs_light", "Heavy W mass vs light W qq", nbins, 0.0, 120.0, nbins, 0.0, 120.0);

    TH1F* met_x_corr = new TH1F("met_x_corr", "MET Px corrected", nbins, 0.0, 350.0);
    TH1F* met_y_corr = new TH1F("met_y_corr", "MET Py corrected", nbins, 0.0, 350.0);
    TH1F* met_x = new TH1F("met_x", "MET Px", nbins, 0.0, 350.0);
    TH1F* met_y = new TH1F("met_y", "MET Py", nbins, 0.0, 350.0);
    THStack* met_px = new THStack("met_px", "MET Px before and after correcting");
    auto met_x_legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    THStack* met_py = new THStack("met_py", "MET Py before and after correcting");
    auto met_y_legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    TH1F* w_from_lv_comp = new TH1F("w_from_lv_comp", "W mass from gen lepton and reconstructed nu", nbins, 0.0, 350.0);

    int nEvents = myTree->GetEntries();

    TwoObjFunc inv_mass = [](TLorentzVector const& p1, TLorentzVector const& p2) -> float { return (p1 + p2).M(); };
    TwoObjFunc pdf = [](TLorentzVector const& q, TLorentzVector const& j) -> float { return q.Pt()/j.Pt(); };

    int n_onshell = 0;
    int anom = 0;

    TRandom3 rg;

    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector met;
        met.SetPtEtaPhiM(genMET_pT, 0, genMET_phi, 0);

        met_x->Fill(met.Px());
        met_y->Fill(met.Py());

        TLorentzVector q1, q2, j1, j2, l, nu;
        q1.SetPtEtaPhiM(genq1_pt, genq1_eta, genq1_phi, genq1_mass);
        q2.SetPtEtaPhiM(genq2_pt, genq2_eta, genq2_phi, genq2_mass);
        j1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        j2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);

        genW_vs_lv->Fill(inv_mass(l, nu), genw_mass);
        genW_vs_qq->Fill(inv_mass(q1, q2), genw_mass);

        fill_2obj_dist(l, nu, inv_mass, w_from_lv);
        h_from_qqlv->Fill((q1+q2+l+nu).M());
        // if ((q1+q2+l+nu).M() > 130.0)
        // {
        //     std::cout << "q1: ";
        //     q1.Print();
        //     std::cout << "q2: ";
        //     q2.Print();
        //     std::cout << "l: ";
        //     l.Print();
        //     std::cout << "nu: ";
        //     nu.Print();
        //     std::cout << "----------------------------------------------------------" << std::endl;
        // }

        if ((l+nu).M() > (q1+q2).M())
        {
            // offshell_w_from_qq->Fill((q1+q2).M());
            // if ((q1+q2).M() > 75)
            // {
            //     std::cout << "m_qq = " << (q1+q2).M() << ", m_lv = " 
            //                            << (l+nu).M() << ", sum = " 
            //                            << (l+nu).M() + (q1+q2).M() << ", m_h = " 
            //                            << (q1+q2+l+nu).M() << std::endl;
            //     ++anom;
            // }
            fill_2obj_dist(q1, q2, inv_mass, offshell_w_from_qq);
            heavy_vs_light->Fill(inv_mass(l, nu), inv_mass(q1, q2));
            // std::cout << inv_mass(q1, q2) << std::endl;
            // fill_2obj_dist(l, nu, inv_mass, onshell_w_from_qq);
        }
        else
        {
            // onshell_w_from_qq->Fill((q1+q2).M());
            fill_2obj_dist(q1, q2, inv_mass, onshell_w_from_qq);
            heavy_vs_light->Fill(inv_mass(q1, q2), inv_mass(l, nu));
            // fill_2obj_dist(l, nu, inv_mass, offshell_w_from_qq);
        }

        if ((l+nu).M() > (j1+j2).M())
        {
            fill_2obj_dist(j1, j2, inv_mass, offshell_w_from_jj);
            // fill_2obj_dist(l, nu, inv_mass, onshell_w_from_jj);
        }
        else
        {
            fill_2obj_dist(j1, j2, inv_mass, onshell_w_from_jj);
            // fill_2obj_dist(l, nu, inv_mass, offshell_w_from_jj);
        }

        pt_order(q1, q2);
        pt_order(j1, j2);
        fill_2obj_dist(q1, q2, inv_mass, w_from_qq);
        fill_2obj_dist(j1, j2, inv_mass, w_from_jj);

        // if (is_offshell(j1, j2, W_threshold))
        // {
        //     fill_2obj_dist(q1, j1, pdf, lead_off);
        // }
        // else
        // {
        //     fill_2obj_dist(q1, j1, pdf, lead_on);
        // }
        if (is_offshell(q1, q2, l, nu))
        {
            fill_2obj_dist(q1, j1, pdf, lead_off);
        }
        else
        {
            fill_2obj_dist(q1, j1, pdf, lead_on);
        }

        TLorentzVector bj1, bj2, bq1, bq2;
        bj1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        bj2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);
        pt_order(bj1, bj2);
        pt_order(bq1, bq2);
        fill_2obj_dist(bq1, bj1, pdf, lead_bjet_pdf);
        // fill_2obj_dist(bq1, bq2, inv_mass, h_from_bq);
        // TwoObjFunc lambda = [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaR(p2); };
        // fill_2obj_dist(bq1, bj1, lambda, bjbq_dR);
        // lead_q_pt->Fill(bq1.Pt());
        // lead_j_pt->Fill(bj1.Pt());
    }

    // std::cout << "# of cases with (q1+q2+l+nu).M() > m_h: " << anom << std::endl;

    for (size_t i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector b1, b2, p1, p2, l, nu, met;
        b1.SetPtEtaPhiM(genbjet1_pt, genbjet1_eta, genbjet1_phi, genbjet1_mass);
        b2.SetPtEtaPhiM(genbjet2_pt, genbjet2_eta, genbjet2_phi, genbjet2_mass);
        l.SetPtEtaPhiM(genl_pt, genl_eta, genl_phi, genl_mass);
        nu.SetPtEtaPhiM(gennu_pt, gennu_eta, gennu_phi, gennu_mass);
        met.SetPtEtaPhiM(genMET_pT, rg.Uniform(-6, 6), genMET_phi, 0);

        // fill_2obj_dist(p1, p2, inv_mass, h_from_bb);

        float mh = 125.0;
        std::pair<float, float> b_resc_factors = compute_jet_resc_factors(b1, b2, lead_bjet_pdf, mh);
        // rescale_jets(p1, p2, lead_bjet_pdf, mh);
        // rescale_jets(p1, p2, b_jet_rescale, rg, mh);
        // lead_j_pt_corr->Fill(p1.Pt());
        // fill_2obj_dist(p1, p2, inv_mass, h_from_bb_corr);

        p1.SetPtEtaPhiM(genqjet1_pt, genqjet1_eta, genqjet1_phi, genqjet1_mass);
        p2.SetPtEtaPhiM(genqjet2_pt, genqjet2_eta, genqjet2_phi, genqjet2_mass);

        pt_order(p1, p2);
        pt_order(b1, b2);

        // if (is_offshell(p1, p2, W_threshold))
        std::pair<float, float> light_j_resc_factors;
        if (is_offshell(p1, p2, l, nu))
        {
            // rescale_jets(p1, p2, lead_off, off_shell_mass);
            // rescale_jets(p1, p2, lead_off, offshell_w_from_qq);
            light_j_resc_factors = compute_jet_resc_factors(p1, p2, lead_off, offshell_w_from_qq);
        }
        else
        {
            // rescale_jets(p1, p2, lead_on, on_shell_mass);
            // rescale_jets(p1, p2, lead_on, onshell_w_from_qq);
            light_j_resc_factors = compute_jet_resc_factors(p1, p2, lead_on, onshell_w_from_qq);
        }
        TLorentzVector met_corr = compute_met_corr(met, b1, b2, p1, p2, b_resc_factors, light_j_resc_factors, rg, 25.2);
        met_x_corr->Fill(met_corr.Px());
        met_y_corr->Fill(met_corr.Py());

        TLorentzVector nu_comp = generate_nu(met_corr);
        fill_2obj_dist(nu_comp, l, inv_mass, w_from_lv_comp);

        p1 *= light_j_resc_factors.first;
        p2 *= light_j_resc_factors.second;
        fill_2obj_dist(p1, p2, inv_mass, w_from_jj_corr);
    }
    // std::cout << "Onshell W bosons: " << n_onshell << std::endl;

    // std::cout << "b jet corrections failed: " << fails << std::endl;
    // std::cout << "b jet corrections total: " << total << std::endl;
    // std::cout << "b jet correction effieciency: " <<  "1 - " <<  fails << "/" << total << " = " << 1 - fails*1.0/total << std::endl;

    // std::cout << "q jet corrections failed: " << fails << std::endl;
    // std::cout << "q jet corrections total: " << total << std::endl;
    // std::cout << "q jet correction effieciency: " <<  "1 - " <<  fails << "/" << total << " = " << 1 - fails*1.0/total << std::endl;
    
    // save_1d_dist(w_from_qq, "w_from_qq", "[GeV]");
    // save_1d_dist(h_from_bb_corr, "h_from_bb_corr", "[GeV]");
    // save_1d_dist(h_from_bb, "h_from_bb", "[GeV]");
    // save_1d_dist(lead_bjet_pdf, "lead_bjet_pdf", "quark_pt/jet_pt");
    save_1d_dist(lead_on, "lead_on", "quark_pt/jet_pt");
    save_1d_dist(lead_off, "lead_off", "quark_pt/jet_pt");
    save_1d_dist(w_from_jj_corr, "w_from_jj_corr_wpdf", "[GeV]");

    std::vector<TH1F*> hists = {w_from_qq, w_from_jj, w_from_jj_corr};
    std::vector<std::string> legends = {"W mass from qq", "W mass from jj", "W mass from corr jj"};
    save_1d_stack(hists, legends, "corrections", "W mass", "[GeV]");

    save_2d_dist(genW_vs_lv, "genW_vs_lv", "lv mass [GeV]", "gen W mass [GeV]");
    save_2d_dist(genW_vs_qq, "genW_vs_qq", "qq mass [GeV]", "gen W mass [GeV]");

    save_2d_dist(heavy_vs_light, "heavy_vs_light", "heavier W mass [GeV]", "lighter W mass [GeV]");

    hists = {met_x, met_x_corr};
    legends = {"MET Px", "MET Px corrected"};
    save_1d_stack(hists, legends, "met_x_corrections", "MET Px corrections", "[GeV]");

    hists = {met_y, met_y_corr};
    legends = {"MET Py", "MET Py corrected"};
    save_1d_stack(hists, legends, "met_y_corrections", "MET Py corrections", "[GeV]");

    save_1d_dist(w_from_lv_comp, "w_from_lv_comp", "[GeV]");

    // save_1d_dist(lead_j_pt, "lead_j_pt", "[GeV]");
    // save_1d_dist(lead_j_pt_corr, "lead_j_pt_corr", "[GeV]");
    // save_1d_dist(lead_q_pt, "lead_q_pt", "[GeV]");

    // save_1d_dist(bjbq_dR, "bjbq_dR", "dR");

    // save_1d_dist(h_from_bq, "h_from_bq", "[GeV]");
    // save_1d_dist(w_from_lv, "w_from_lv", "[GeV]");

    // save_1d_dist(h_from_qqlv, "h_from_qqlv", "[GeV]");

    // save_1d_dist(onshell_w_from_qq, "onshell_w_from_qq", "[GeV]");
    // save_1d_dist(onshell_w_from_jj, "onshell_w_from_jj", "[GeV]");
    // save_1d_dist(offshell_w_from_qq, "offshell_w_from_qq", "[GeV]");
    // save_1d_dist(offshell_w_from_jj, "offshell_w_from_jj", "[GeV]");

    // std::vector<TH1F*> hists = {lead_j_pt, lead_j_pt_corr, lead_q_pt};
    // std::vector<std::string> legends = {"jet pt", "jet pt corrected", "quark pt"};
    // save_1d_stack(hists, legends, "pts", "Leading b object Pt", "[GeV]");

    return 0;
}