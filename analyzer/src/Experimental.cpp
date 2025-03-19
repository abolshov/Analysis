#include "Experimental.hpp"

namespace Experimental 
{
    namespace SL 
    {
        std::optional<Float_t> ComputeSublRescFromMassConstr(LorentzVectorF_t const& lead_jet,
                                                             LorentzVectorF_t const& subl_jet,
                                                             Float_t lead_resc, Float_t mass)
        {
            Float_t x1 = subl_jet.M2();
            Float_t x2 = 2.0*lead_resc*(lead_jet.Dot(subl_jet));
            Float_t x3 = lead_resc*lead_resc*lead_jet.M2() - mass*mass;
            Float_t discrim = x2*x2 - 4.0*x1*x3;
            if (x2 > 0.0 && x1 != 0.0 && discrim > 0.0)
            {
                Float_t subl_resc = (-x2 + std::sqrt(discrim))/(2.0*x1);
                if (subl_resc > 0.0)
                {
                    return std::make_optional<subl_resc>;
                }
            }
            return std::nullopt;
        }   

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

                LorentzVectorF_t j1 = lj1;
                LorentzVectorF_t j2 = lj2;

                Double_t lead_resc = 1.0;
                Double_t mWhad = 1.0;
                pdf_c3mWhad>GetRandom2(lead_resc, mWhad, prg.get());
                auto subl_resc = ComputeSublRescFromMassConstr(j1, j2, lead_resc, mWhad);
                if (!subl_resc.has_value())
                {
                    continue;
                }

                Float_t c3 = lead_resc;
                Float_t c4 = subl_resc.value();
                Float_t ljet_resc_dpx = -1.0*(c3 - 1)*lj1.Px() - (c4 - 1)*lj2.Px();
                Float_t ljet_resc_dpy = -1.0*(c3 - 1)*lj1.Py() - (c4 - 1)*lj2.Py();
                j1 *= c3;
                j2 *= c4;
            }
        }
    }
}