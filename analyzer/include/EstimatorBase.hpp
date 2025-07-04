#ifndef ESTIM_BASE
#define ESTIM_BASE

#include <optional>
#include <type_traits>

#include "TRandom3.h"

#include "Definitions.hpp"
#include "Recorder.hpp"
#include "Constants.hpp"
#include "JetCombination.hpp"
#include "Event.hpp"

template <typename T, std::enable_if_t<std::is_default_constructible_v<T>, bool> = true>
void ResetObject(T& object)
{
    object = T{};
}

class EstimatorBase
{
    public:
    explicit EstimatorBase(AggregationMode aggr_mode);
    EstimatorBase(TString const& dbg_file_name, AggregationMode aggr_mode);
    virtual ~EstimatorBase() = default;

    // EstimateMass must implement logic of estimation: it should decide what data to use for estimation based on what's available 
    virtual OptArrF_t<ESTIM_OUT_SZ> EstimateMass(Event const& event) = 0;
    inline void OpenDbgFile(TString const& dbg_file_name) { m_recorder.OpenFile(dbg_file_name); }

    protected:
    HistVec_t<TH1F> m_pdf_1d;
    HistVec_t<TH2F> m_pdf_2d;
    std::unique_ptr<TRandom3> m_prg;
    UHist_t<TH1F> m_res_mass;
    Recorder m_recorder; 
    AggregationMode m_aggr_mode;

    virtual std::unique_ptr<TTree> MakeTree(TString const& tree_name) = 0;
};

#endif