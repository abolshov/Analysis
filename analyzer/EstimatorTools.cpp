#include "EstimatorTools.hpp"

LorentzVectorF_t SamplePNetResCorr(LorentzVectorF_t const& jet, std::unique_ptr<TRandom3>& prg, Float_t resolution)
{
    Float_t dpt = prg->Gaus(0.0, resolution);
    Float_t pt = jet.Pt();
    while (pt + dpt < 0.0)
    {
        dpt = prg->Gaus(0.0, resolution);
    }
    return LorentzVectorF_t(pt + dpt, jet.Eta(), jet.Phi(), jet.M());
}

std::pair<Float_t, Float_t> ComputeJetResc(LorentzVectorF_t const& b1, LorentzVectorF_t const& b2, UHist_t<TH1F>& pdf, Float_t mass)
{
    Float_t c1 = pdf->GetRandom();
    Float_t x1 = b2.M2();
    Float_t x2 = 2.0*c1*(b1.Dot(b2));
    Float_t x3 = c1*c1*b1.M2() - mass*mass;
    Float_t discrim = x2*x2 - 4.0*x1*x3;
    Float_t c2 = 1.0;
    if (x2 >= 0.0 && x1 != 0.0 && discrim >= 0.0)
    {
        c2 = (-x2 + std::sqrt(x2*x2 - 4.0*x1*x3))/(2.0*x1);
        if (c2 >= 0.0)
        {
            return std::make_pair(c1, c2); 
        }
    }
    return {1.0, 1.0};
}

std::optional<LorentzVectorF_t> NuFromOnshellW(Float_t eta, Float_t phi, Float_t mw, LorentzVectorF_t const& lep_onshell)
{
    Float_t deta = eta - lep_onshell.Eta();
    Float_t dphi = phi - lep_onshell.Phi();
    Float_t pt = mw*mw/(2.0*lep_onshell.Pt()*(std::cosh(deta) - std::cos(dphi)));

    if (std::isinf(pt) || std::isnan(pt))
    {
        return std::nullopt;
    }
    return std::make_optional<LorentzVectorF_t>(pt, eta, phi, 0.0);
}

std::optional<LorentzVectorF_t> NuFromOffshellW(LorentzVectorF_t const& lep1, 
                                                LorentzVectorF_t const& lep2, 
                                                LorentzVectorF_t const& nu1,
                                                LorentzVectorF_t const& met,
                                                int control, 
                                                Float_t mh)
{
    LorentzVectorF_t tmp = lep1 + lep2 + nu1;

    // how much is contribution of nu2 in met
    // assumption: met = nu1 + nu2
    Float_t nu_tmp_px = met.Px() - nu1.Px();
    Float_t nu_tmp_py = met.Py() - nu1.Py();

    // weird stuff 
    Float_t px = std::sqrt(tmp.Pt()*tmp.Pt() + tmp.M2());
    Float_t py = 0.0;
    Float_t pz = tmp.Pz();
    Float_t en = tmp.E();

    LorentzVectorF_t tmp2;
    tmp2.SetPxPyPzE(px, py, pz, en);

    // pt and phi of nu2 are determined by difference of met and transverse component of nu1
    Float_t pt = std::sqrt(nu_tmp_px*nu_tmp_px + nu_tmp_py*nu_tmp_py);
    Float_t phi = std::atan2(nu_tmp_px, nu_tmp_py);

    Float_t cosh_deta = (mh*mh + nu_tmp_px*tmp.Px() + nu_tmp_py*tmp.Py() - tmp.M2())/(2.0*tmp2.Pt()*pt);
    if (cosh_deta < 1.0)
    {
        return std::nullopt;
    }

    Float_t delta_eta = std::acosh(cosh_deta);
    Float_t eta = control == 1 ? tmp2.Eta() - delta_eta : tmp2.Eta() + delta_eta;

    if (std::abs(eta) > 7.0)
    {
        return std::nullopt;
    }

    return std::make_optional<LorentzVectorF_t>(pt, eta, phi, 0.0);
}