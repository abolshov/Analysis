#ifndef EVAL_DL_HPP
#define EVAL_DL_HPP

#include "EvaluatorBase.hpp"

class EvaluatorDoubleLep : public EvaluatorBase
{
    public:
    using EvaluatorBase::EvaluatorBase;
    ~EvaluatorDoubleLep() = default;
    
    void Evaluate() override;
};

#endif