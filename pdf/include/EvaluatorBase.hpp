#ifndef EVAL_BASE_HPP
#define EVAL_BASE_HPP

#include <memory>
#include <fstream>

#include "Definitions.hpp"
#include "Buffer.hpp"

#include "TTree.h"
#include "TChain.h"

class EvaluatorBase
{
    public:
    EvaluatorBase(std::ifstream& file_stream, Channel ch);
    virtual ~EvaluatorBase() = default;

    virtual void Evaluate() = 0;

    protected:
    std::unique_ptr<TChain> m_chain;
    std::unique_ptr<Buffer> m_buf;
};

#endif