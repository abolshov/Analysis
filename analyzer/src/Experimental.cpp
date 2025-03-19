#include "Experimental.hpp"

namespace Experimental 
{
    namespace SL 
    {
        ArrF_t<ESTIM_OUT_SZ> EstimateCombination(VecLVF_t const& particles, 
                                                 HistVec_t<TH1F> const& pdfs_1d,
                                                 HistVec_t<TH2F> const& pdfs_2d, 
                                                 std::unique_ptr<TTree>& dbg_tree,
                                                 std::unique_ptr<TRandom3>& prg,
                                                 ULong64_t evt_id, 
                                                 TString const& comb_label)
        {
            ArrF_t<ESTIM_OUT_SZ> res{};
            std::fill(res.begin(), res.end(), -1.0f);

            LorentzVectorF_t const& bj1 = particles[static_cast<size_t>(ObjSL::bj1)];
            LorentzVectorF_t const& bj2 = particles[static_cast<size_t>(ObjSL::bj2)];
            LorentzVectorF_t const& lj1 = particles[static_cast<size_t>(ObjSL::lj1)];
            LorentzVectorF_t const& lj2 = particles[static_cast<size_t>(ObjSL::lj2)];
            LorentzVectorF_t const& lep = particles[static_cast<size_t>(ObjSL::lep)];
            LorentzVectorF_t const& met = particles[static_cast<size_t>(ObjSL::met)];

            UHist_t<TH2F>& pdf_c1c2 = m_pdf_2d[static_cast<size_t>(PDF2_sl::c1c2)];
            UHist_t<TH2F>& pdf_c3mWhad = m_pdf_2d[static_cast<size_t>(PDF2_sl::c3mWhad)];

            TString title = Form("X->HH mass: event %llu, comb %s", evt_id, comb_label.Data());
            auto res_mass = std::make_unique<TH1F>("X_mass", title, N_BINS, MIN_MASS, MAX_MASS);

            Float_t mh = prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);

            for (int i = 0; i < N_ITER; ++i)
            {
                Float_t smear_dpx = prg->Gaus(0.0, MET_SIGMA);
                Float_t smear_dpy = prg->Gaus(0.0, MET_SIGMA);

                Double_t c1 = 1.0;
                Double_t c2 = 1.0;
                pdf_c1c2->GetRandom2(c1, c2, prg.get());

                LorentzVectorF_t b1 = bj1;
                LorentzVectorF_t b2 = bj2;
                b1 *= c1;
                b2 *= c2;

                LorentzVectorF_t Hbb = b1 + b2;
                if (std::abs(Hbb.M() - mh) > 1.0)
                {
                    continue;
                }

                Float_t bjet_resc_dpx = -1.0*(c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px();
                Float_t bjet_resc_dpy = -1.0*(c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py();
            }
        }
    }
}