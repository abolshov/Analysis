#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>

#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
// #include "TH1D.h"
// #include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

using RescFact = std::pair<float, float>;
struct CorrData 
{
    RescFact resc;
    TLorentzVector p1;
    TLorentzVector p2;
};

void Print(TLorentzVector const& p, bool EXYZ = false);

namespace save 
{
    void save_1d_dist(TH1F* dist, std::string const& name, std::string const& title);

    void save_2d_dist(TH2F* dist, std::string const& name, std::string const& title_x, std::string const& title_y);

    void save_1d_stack(std::vector<TH1F*> const& distrs,
                       std::vector<std::string> const& legends,
                       std::string const& name,
                       std::string const& title,
                       std::string const& axis_label); 

    std::pair<double, double> save_fit(TH1F* dist, std::string const& name, std::string const& title, std::string const& fit_func);   
}

namespace jet
{
    void pt_order(TLorentzVector& p1, TLorentzVector& p2);

    bool is_offshell(TLorentzVector const& p1, TLorentzVector const& p2, TLorentzVector const& l, TLorentzVector const& nu);

    std::pair<float, float> compute_resc_factors(TLorentzVector& p1, TLorentzVector& p2, TH1F* lead_jet_pdf, TH1F* mass_pdf);
}

TLorentzVector NuFromLeptonicW_v1(float nu_eta, float nu_phi, TLorentzVector const& l, float mw);
TLorentzVector NuFromLeptonicW_v2(TLorentzVector const& l, TLorentzVector const& j1, TLorentzVector const& j2, TLorentzVector const& met, float mh);

namespace hme
{
    enum GEN_PART { b1, b2, j1, j2, l, nu, met };
    // enum PDF {};

    float analytical(std::vector<TLorentzVector> const& particles);
    float rand_sampl(std::vector<TLorentzVector> const& particles, std::vector<TH1F*> const& pdfs, int nIter, TRandom3& rg, int nbins, bool correct_light_jets = false, bool weighted = false);
}

namespace kinematics
{
    float mT(TLorentzVector const& p1, TLorentzVector const& p2);
}

bool HasZeroParticle(std::vector<TLorentzVector> const& particles);
bool HasIdenticalPair(TLorentzVector const& p1, TLorentzVector const& p2);
bool ValidDeltaR(TLorentzVector const& p1, TLorentzVector const& p2);

#endif