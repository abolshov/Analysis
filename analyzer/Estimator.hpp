#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

#include "TRandom3.h"

#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "Constants.hpp"

class EstimatorSingLep
{
    public:
    EstimatorSingLep() = default;
    explicit EstimatorSingLep(TString const& file_name);

    std::array<Float_t, OUTPUT_SIZE> EstimateMass(VecLVF_t const& particles, ULong64_t evt, int comb_id, std::pair<Float_t, Float_t> lj_pt_res);

    private:
    HistVec_t<TH1F> pdf_1d;
    HistVec_t<TH2F> pdf_2d;
    std::unique_ptr<TRandom3> prg;
};

#endif