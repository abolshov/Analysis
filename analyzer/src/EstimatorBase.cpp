#include "EstimatorBase.hpp"

EstimatorBase::EstimatorBase() 
:   m_prg(std::make_unique<TRandom3>(SEED))
,   m_res_mass(std::make_unique<TH1F>("none", "none", N_BINS, MIN_MASS, MAX_MASS))
{}

EstimatorBase::EstimatorBase(TString const& dbg_file_name) 
:   m_prg(std::make_unique<TRandom3>(SEED))
,   m_res_mass(std::make_unique<TH1F>("none", "none", N_BINS, MIN_MASS, MAX_MASS))
,   m_recorder(dbg_file_name)
{}
