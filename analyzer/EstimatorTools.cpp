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

std::pair<Float_t, Float_t> ComputeJetResc(LorentzVectorF_t const& b1, LorentzVectorF_t const& b2, UHist_t& pdf, Float_t mass)
{
    Float_t c1 = pdf->GetRandom();
    Float_t x1 = b2.M2();
    Float_t x2 = 2.0*c1*(b1.Dot(b2));
    Float_t x3 = c1*c1*b1.M2() - mh*mh;
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
