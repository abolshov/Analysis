#ifndef EXP_HPP
#define EXP_HPP

#include <optional>

#include "Constants.hpp"
#include "Definitions.hpp"
#include "EstimatorBase.hpp"
#include "EstimatorTools.hpp"
#include "EstimatorUtils.hpp"

namespace Experimental 
{
    namespace SL
    {
        std::optional<Float_t> ComputeSublRescFromMassConstr(LorentzVectorF_t const& lead_jet,
                                                             LorentzVectorF_t const& subl_jet,
                                                             Float_t lead_resc, Float_t mass);

        Float_t GenW2Mass(UHist_t<TH1F>& pdf, std::unique_ptr<TRandom3>& prg, Float_t mh, Float_t mw1);
        
        ArrF_t<ESTIM_OUT_SZ> EstimateCombination(VecLVF_t const& particles, 
                                                 HistVec_t<TH1F>& pdfs_1d,
                                                 HistVec_t<TH2F>& pdfs_2d, 
                                                 UHist_t<TH1F>& res_mass,
                                                 std::unique_ptr<TRandom3>& prg,
                                                 ULong64_t evt_id, 
                                                 TString const& comb_label);
    }
}

#endif