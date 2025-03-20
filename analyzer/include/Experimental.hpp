#ifndef EXP_HPP
#define EXP_HPP

#include <optional>

#include "Constants.hpp"
#include "Definitions.hpp"
#include "EstimatorBase.hpp"

namespace Experimental 
{
    namespace SL
    {
        std::optional<Float_t> ComputeSublRescFromMassConstr(LorentzVectorF_t const& lead_jet,
                                                             LorentzVectorF_t const& subl_jet,
                                                             Float_t lead_resc, Float_t mass);

        Float_t GenLepWMass(UHist_t<TH1F>& pdf, std::unique_ptr<TRandom3>& prg, Float_t mh, Float_t mWhad);
        
        ArrF_t<ESTIM_OUT_SZ> EstimateCombination(VecLVF_t const& particles, 
                                                 HistVec_t<TH1F> const& pdfs_1d,
                                                 HistVec_t<TH2F> const& pdfs_2d, 
                                                 std::unique_ptr<TTree>& dbg_tree,
                                                 std::unique_ptr<TRandom3>& prg,
                                                 ULong64_t evt_id, 
                                                 TString const& comb_label);
    }
}

#endif