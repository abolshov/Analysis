#include "EvaluatorDoubleLep.hpp"
#include "Definitions.hpp"
#include "Constants.hpp"

#include <iostream>
#include <memory>

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

#include "TH1F.h"

bool CorrectLepReco(int lep_type, int lep_genLep_kind)
{
    bool reco_lep_mu = (lep_type == 2);
    bool reco_lep_ele = (lep_type == 1);

    bool gen_lep_mu = ((lep_genLep_kind == 2) || (lep_genLep_kind == 4));
    bool gen_lep_ele = ((lep_genLep_kind == 1) || (lep_genLep_kind == 3));

    bool corr_lep_reco = ((reco_lep_mu && gen_lep_mu) || (reco_lep_ele && gen_lep_ele));
    return corr_lep_reco;
}

void EvaluatorDoubleLep::Evaluate()
{
    auto pdf_fatbb = std::make_unique<TH1F>("pdf_fatbb", "1d PDF for fat H->bb correction", N_BINS, 0, 8);
    ULong64_t n_events = m_chain->GetEntries();
    std::cout << "Number of events: " << n_events << "\n";
    for (ULong64_t i = 0; i < n_events; ++i)
    {
        m_chain->GetEntry(i);

        if (!(m_buf->reco_lep_pt[static_cast<size_t>(Lep::lep1)] > 0.0f && m_buf->reco_lep_pt[static_cast<size_t>(Lep::lep2)] > 0.0f))
        {
            continue;
        }

        bool correct_lep1 = CorrectLepReco(m_buf->reco_lep_type[static_cast<size_t>(Lep::lep1)], m_buf->reco_lep_gen_kind[static_cast<size_t>(Lep::lep1)]);
        bool correct_lep2 = CorrectLepReco(m_buf->reco_lep_type[static_cast<size_t>(Lep::lep2)], m_buf->reco_lep_gen_kind[static_cast<size_t>(Lep::lep2)]);
        bool correct_leptons = correct_lep1 && correct_lep2;
        if (!correct_leptons)
        {
            continue;
        }

        if (m_buf->n_reco_fatjet < 1)
        {
            continue;
        }

        if (m_buf->gen_quark_pt[static_cast<size_t>(Quark::b1)] < 20.0f || m_buf->gen_quark_pt[static_cast<size_t>(Quark::b2)] < 20.0f)
        {
            continue;
        }

        if (std::abs(m_buf->gen_quark_eta[static_cast<size_t>(Quark::b1)]) > 2.5f || std::abs(m_buf->gen_quark_eta[static_cast<size_t>(Quark::b2)]) > 2.5f)
        {
            continue;
        }
        
        LorentzVectorF_t b1_p4(m_buf->gen_quark_pt[static_cast<size_t>(Quark::b1)], 
                               m_buf->gen_quark_eta[static_cast<size_t>(Quark::b1)],
                               m_buf->gen_quark_phi[static_cast<size_t>(Quark::b1)],
                               m_buf->gen_quark_mass[static_cast<size_t>(Quark::b1)]);
        LorentzVectorF_t b2_p4(m_buf->gen_quark_pt[static_cast<size_t>(Quark::b2)], 
                               m_buf->gen_quark_eta[static_cast<size_t>(Quark::b2)],
                               m_buf->gen_quark_phi[static_cast<size_t>(Quark::b2)],
                               m_buf->gen_quark_mass[static_cast<size_t>(Quark::b2)]);
        if (DeltaR(b1_p4, b2_p4) > 0.4)
        {
            continue;
        }
        
        LorentzVectorF_t hbb_p4(m_buf->hbb_pt, m_buf->hbb_eta, m_buf->hbb_phi, m_buf->hbb_mass);
        LorentzVectorF_t reco_hbb_p4{};
        size_t num_fatjet = m_buf->n_reco_fatjet;
        bool have_fatbb_match = false;
        Float_t min_dr = 10.0f;
        for (size_t i = 0; i < num_fatjet; ++i)
        {
            LorentzVectorF_t tmp = LorentzVectorF_t(m_buf->reco_fatjet_pt[i], 
                                                    m_buf->reco_fatjet_eta[i],
                                                    m_buf->reco_fatjet_phi[i],
                                                    m_buf->reco_fatjet_mass[i]);
            Float_t dr = DeltaR(tmp, hbb_p4);
            if (dr < 0.4 && dr < min_dr)
            {
                reco_hbb_p4 = tmp;
                have_fatbb_match = true;
            }
        }

        if (!have_fatbb_match)
        {
            continue;
        }

        pdf_fatbb->Fill(m_buf->hbb_pt/reco_hbb_p4.Pt());
    }

    int binmax = pdf_fatbb->GetMaximumBin();
    Float_t peak_val = pdf_fatbb->GetBinContent(binmax);
    pdf_fatbb->Scale(1.0/peak_val);
    pdf_fatbb->Write();
	m_out_file->Write();
	m_out_file->Close();
}