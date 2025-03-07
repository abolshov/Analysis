#include "EstimatorUtils.hpp"

VecLVF_t GetRecoJetP4(Storage const& s)
{
    VecLVF_t res;
    for (int i = 0; i < s.n_reco_jet; ++i)
    {
        res.emplace_back(s.reco_jet_pt[i], s.reco_jet_eta[i], s.reco_jet_phi[i], s.reco_jet_mass[i]);
    }
    return res;
}

VecLVF_t GetRecoLepP4(Storage const& s, Channel ch)
{
    VecLVF_t res;
    res.emplace_back(s.reco_lep_pt[static_cast<size_t>(Lep::lep1)], 
                     s.reco_lep_eta[static_cast<size_t>(Lep::lep1)], 
                     s.reco_lep_phi[static_cast<size_t>(Lep::lep1)], 
                     s.reco_lep_mass[static_cast<size_t>(Lep::lep1)]);
    if (ch == Channel::DL)
    {
        res.emplace_back(s.reco_lep_pt[static_cast<size_t>(Lep::lep2)], 
                         s.reco_lep_eta[static_cast<size_t>(Lep::lep2)], 
                         s.reco_lep_phi[static_cast<size_t>(Lep::lep2)], 
                         s.reco_lep_mass[static_cast<size_t>(Lep::lep2)]);
    }
    return res;
}

std::vector<Float_t> GetPNetRes(Storage const& s)
{
    std::vector<Float_t> res;
    for (int i = 0; i < s.n_reco_jet; ++i)
    {   
        Float_t mult = s.reco_jet_corr[i]*s.reco_jet_res[i] == 0 ? DEFAULT_JET_RES : s.reco_jet_corr[i]*s.reco_jet_res[i];
        res.push_back(s.reco_jet_pt[i]*mult);
    }
    return res;
}

void Get1dPDFs(TFile* fptr, HistVec_t<TH1F>& pdfs, Channel ch)
{
    if (ch == Channel::SL)
    {
        for (auto const& [pdf, name]: pdf1d_sl_names)
        {
            pdfs[static_cast<size_t>(pdf)] = ReadHist<TH1F>(fptr, name);
        }
    }
    else if (ch == Channel::DL)
    {
        for (auto const& [pdf, name]: pdf1d_dl_names)
        {
            pdfs[static_cast<size_t>(pdf)] = ReadHist<TH1F>(fptr, name);
        }
    }
    else 
    {
        throw std::runtime_error("Get1dPDFs: attempting to read PDFs for unnkown channel");
    }
}

void Get2dPDFs(TFile* fptr, HistVec_t<TH2F>& pdfs, Channel ch)
{
    if (ch == Channel::SL)
    {
        for (auto const& [pdf, name]: pdf2d_sl_names)
        {
            pdfs[static_cast<size_t>(pdf)] = ReadHist<TH2F>(fptr, name);
        }
    }
    else if (ch == Channel::DL)
    {
        for (auto const& [pdf, name]: pdf2d_dl_names)
        {
            pdfs[static_cast<size_t>(pdf)] = ReadHist<TH2F>(fptr, name);
        }
    }
    else 
    {
        throw std::runtime_error("Get2dPDFs: attempting to read PDFs for unnkown channel");
    }
}

Float_t ComputeWidth(UHist_t<TH1F> const& h, unsigned l, unsigned r)
{
    int const nq = 100;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    Double_t yq[nq];  // array to contain the quantiles
    for (int i = 0; i < nq; ++i) 
    {
        xq[i] = static_cast<Double_t>(i+1)/nq;
    }
    h->GetQuantiles(nq, yq, xq);
    return static_cast<Float_t>(yq[r - 1] - yq[l - 1]);
}

#ifdef DEBUG
VecLVF_t GetGenJetP4(Storage const& s)
{
    VecLVF_t res;
    for (int i = 0; i < s.n_gen_jet; ++i)
    {
        res.emplace_back(s.gen_jet_pt[i], s.gen_jet_eta[i], s.gen_jet_phi[i], s.gen_jet_mass[i]);
    }
    return res;
}

VecLVF_t GetGenLepP4(Storage const& s, Channel ch)
{
    VecLVF_t res;
    res.emplace_back(s.gen_lep_pt[static_cast<size_t>(Lep::lep1)], 
                     s.gen_lep_eta[static_cast<size_t>(Lep::lep1)], 
                     s.gen_lep_phi[static_cast<size_t>(Lep::lep1)], 
                     s.gen_lep_mass[static_cast<size_t>(Lep::lep1)]);
    if (ch == Channel::DL)
    {
        res.emplace_back(s.gen_lep_pt[static_cast<size_t>(Lep::lep2)], 
                         s.gen_lep_eta[static_cast<size_t>(Lep::lep2)], 
                         s.gen_lep_phi[static_cast<size_t>(Lep::lep2)], 
                         s.gen_lep_mass[static_cast<size_t>(Lep::lep2)]);
    }
    return res;
}

VecLVF_t GetGenQuarksP4(Storage const& s, Channel ch)
{
    VecLVF_t res;
    res.emplace_back(s.gen_quark_pt[static_cast<size_t>(Quark::b1)], 
                     s.gen_quark_eta[static_cast<size_t>(Quark::b1)], 
                     s.gen_quark_phi[static_cast<size_t>(Quark::b1)], 
                     s.gen_quark_mass[static_cast<size_t>(Quark::b1)]);
    res.emplace_back(s.gen_quark_pt[static_cast<size_t>(Quark::b2)], 
                     s.gen_quark_eta[static_cast<size_t>(Quark::b2)], 
                     s.gen_quark_phi[static_cast<size_t>(Quark::b2)], 
                     s.gen_quark_mass[static_cast<size_t>(Quark::b2)]);
    if (ch == Channel::SL)
    {
        res.emplace_back(s.gen_quark_pt[static_cast<size_t>(Quark::q1)], 
                         s.gen_quark_eta[static_cast<size_t>(Quark::q1)], 
                         s.gen_quark_phi[static_cast<size_t>(Quark::q1)], 
                         s.gen_quark_mass[static_cast<size_t>(Quark::q1)]);
        res.emplace_back(s.gen_quark_pt[static_cast<size_t>(Quark::q2)], 
                         s.gen_quark_eta[static_cast<size_t>(Quark::q2)], 
                         s.gen_quark_phi[static_cast<size_t>(Quark::q2)], 
                         s.gen_quark_mass[static_cast<size_t>(Quark::q2)]);
    }
    return res;
}

VecLVF_t GetGenNuP4(Storage const& s, Channel ch)
{
    VecLVF_t res;
    res.emplace_back(s.gen_nu_pt[static_cast<size_t>(Nu::nu1)], 
                     s.gen_nu_eta[static_cast<size_t>(Nu::nu1)], 
                     s.gen_nu_phi[static_cast<size_t>(Nu::nu1)], 
                     s.gen_nu_mass[static_cast<size_t>(Nu::nu1)]);
    if (ch == Channel::DL)
    {
        res.emplace_back(s.gen_nu_pt[static_cast<size_t>(Nu::nu2)], 
                         s.gen_nu_eta[static_cast<size_t>(Nu::nu2)], 
                         s.gen_nu_phi[static_cast<size_t>(Nu::nu2)], 
                         s.gen_nu_mass[static_cast<size_t>(Nu::nu2)]);
    }
    return res;
}

void LogP4(std::stringstream& ss, LorentzVectorF_t const& p4, std::string const& name)
{
    ss << "\t" << name << "=(" << p4.Pt() << ", " << p4.Eta() << ", " << p4.Phi() << ", " << p4.M() << ")\n";
}
#endif