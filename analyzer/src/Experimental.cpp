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
                    return std::make_optional<Float_t>(subl_resc);
                }
            }
            return std::nullopt;
        }  
        
        Float_t GenW2Mass(UHist_t<TH1F>& pdf, std::unique_ptr<TRandom3>& prg, Float_t mh, Float_t mw1)
        {
            Float_t mw2 = pdf->GetRandom(prg.get());
            while (mw2 + mw1 > mh)
            {
                mw2 = pdf->GetRandom(prg.get());
            }
            return mw2;
        }

        ArrF_t<ESTIM_OUT_SZ> EstimateCombination(VecLVF_t const& particles, 
                                                 HistVec_t<TH1F>& pdfs_1d,
                                                 HistVec_t<TH2F>& pdfs_2d, 
                                                 UHist_t<TH1F>& res_mass,
                                                //  std::unique_ptr<TTree>& dbg_tree,
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

            UHist_t<TH2F>& pdf_b1b2 = pdfs_2d[static_cast<size_t>(PDF2_sl::b1b2)];
            UHist_t<TH2F>& pdf_c3mWhad = pdfs_2d[static_cast<size_t>(PDF2_sl::c3mWhad)];

            UHist_t<TH1F>& pdf_mWoffshell = pdfs_1d[static_cast<size_t>(PDF1_sl::mWoffshell)];
            UHist_t<TH1F>& pdf_mWonshell = pdfs_1d[static_cast<size_t>(PDF1_sl::mWonshell)];
            
            TString name = Form("Xhh_mass_%llu_%s", evt_id, comb_label.Data());
            TString title = Form("X->HH mass: event %llu, comb %s", evt_id, comb_label.Data());
            res_mass->SetNameTitle(name, title);
            // auto res_mass = std::make_unique<TH1F>("X_mass", title, N_BINS, MIN_MASS, MAX_MASS);

            Float_t mh = prg->Gaus(HIGGS_MASS, HIGGS_WIDTH);

            for (int i = 0; i < N_ITER; ++i)
            {
                Float_t smear_dpx = prg->Gaus(0.0, MET_SIGMA);
                Float_t smear_dpy = prg->Gaus(0.0, MET_SIGMA);

                Double_t c1 = 1.0;
                Double_t c2 = 1.0;
                pdf_b1b2->GetRandom2(c1, c2, prg.get());

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
                pdf_c3mWhad->GetRandom2(lead_resc, mWhad, prg.get());
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
                
                LorentzVectorF_t hadW = j1 + j2;

                Float_t mWlep_offshell = GenW2Mass(pdf_mWoffshell, prg, mh, mWhad);
                Float_t mWlep_onshell = GenW2Mass(pdf_mWonshell, prg, mh, mWhad);

                Float_t met_corr_px = met.Px() + bjet_resc_dpx + ljet_resc_dpx + smear_dpx;
                Float_t met_corr_py = met.Py() + bjet_resc_dpy + ljet_resc_dpy + smear_dpy;

                Float_t met_corr_pt = std::sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
                Float_t met_corr_phi = std::atan2(met_corr_py, met_corr_px);
                LorentzVectorF_t met_corr = LorentzVectorF_t(met_corr_pt, 0.0, met_corr_phi, 0.0);
                
                std::vector<Float_t> masses;
                std::array<LorentzVectorF_t, CONTROL> lepW{};
                std::array<LorentzVectorF_t, CONTROL> Hww{};
                std::array<LorentzVectorF_t, CONTROL> Xhh{};
                std::array<LorentzVectorF_t, CONTROL> nu{};
                std::array<Bool_t, CONTROL> correct_hww_mass{};
                std::array<Float_t, CONTROL> mass{};
                for (int control = 0; control < CONTROL; ++control)
                {
                    bool lepW_onshell = control / 2;
                    bool add_deta = control % 2;
                    Float_t mWlep = lepW_onshell ? mWlep_offshell : mWlep_onshell;
                    auto opt_nu = NuFromW(lep, met_corr, add_deta, mWlep);
                    if (opt_nu.has_value())
                    {
                        nu[control] = opt_nu.value();
                        lepW[control] = nu[control] + lep;
                        Hww[control] = lepW[control] + hadW;
                        Xhh[control] = Hww[control] + Hbb;
                        mass[control] = Xhh[control].M();
                        correct_hww_mass[control] = (std::abs(mh - Hww[control].M()) < 1.0);
                        if (!correct_hww_mass[control])
                        {
                            continue;
                        }
                        masses.push_back(mass[control]);
                    }
                }

                if (masses.empty())
                {
                    continue;
                }

                Int_t num_sol = masses.size();
                Float_t weight = 1.0/num_sol;
                for (auto mass: masses)
                {
                    res_mass->Fill(mass, weight);
                }
            }

            Float_t integral = res_mass->Integral();
            if (res_mass->GetEntries() && integral > 0.0)
            {
                int binmax = res_mass->GetMaximumBin(); 
                res[static_cast<size_t>(EstimOut::mass)] = res_mass->GetXaxis()->GetBinCenter(binmax);
                res[static_cast<size_t>(EstimOut::peak_value)] = res_mass->GetBinContent(binmax);
                res[static_cast<size_t>(EstimOut::width)] = ComputeWidth(res_mass, Q16, Q84);
                res[static_cast<size_t>(EstimOut::integral)] = integral;
                return res;
            }
            return res;
        }
    }
}