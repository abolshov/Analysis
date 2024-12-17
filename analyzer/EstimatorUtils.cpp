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
        res.push_back(s.reco_jet_pt[i]*s.reco_jet_corr[i]*s.reco_jet_res[i]);
    }
    return res;
}

void Get1dPDFs(TFile* fptr, HistVec_t<TH1F>& pdfs)
{
    pdfs[static_cast<size_t>(PDF1::numet_pt)] = ReadHist<TH1F>(fptr, "pdf_numet_pt");
    pdfs[static_cast<size_t>(PDF1::numet_dphi)] = ReadHist<TH1F>(fptr, "pdf_numet_dphi");
    pdfs[static_cast<size_t>(PDF1::nulep_deta)] = ReadHist<TH1F>(fptr, "pdf_nulep_deta");
    pdfs[static_cast<size_t>(PDF1::hh_dphi)] = ReadHist<TH1F>(fptr, "pdf_hh_dphi");
    pdfs[static_cast<size_t>(PDF1::mbb)] = ReadHist<TH1F>(fptr, "pdf_mbb");
    pdfs[static_cast<size_t>(PDF1::mww)] = ReadHist<TH1F>(fptr, "pdf_mww_narrow");
    pdfs[static_cast<size_t>(PDF1::hh_deta)] = ReadHist<TH1F>(fptr, "pdf_hh_deta");
}

void Get2dPDFs(TFile* fptr, HistVec_t<TH2F>& pdfs)
{
    pdfs[static_cast<size_t>(PDF2::b1b2)] = ReadHist<TH2F>(fptr, "pdf_b1b2");
    pdfs[static_cast<size_t>(PDF2::hh_dEtadPhi)] = ReadHist<TH2F>(fptr, "pdf_hh_dEtadPhi");
    pdfs[static_cast<size_t>(PDF2::hh_pt_e)] = ReadHist<TH2F>(fptr, "pdf_hh_pt_e");
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