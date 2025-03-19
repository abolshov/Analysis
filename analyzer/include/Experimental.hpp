#ifndef EXP_HPP
#define EXP_HPP

#include "Constants.hpp"
#include "Definitions.hpp"
#include "EstimatorBase.hpp"

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
                                                 TString const& comb_label);
    }
}

#endif