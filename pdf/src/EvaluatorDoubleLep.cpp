#include "EvaluatorDoubleLep.hpp"

#include <iostream>

void EvaluatorDoubleLep::Evaluate()
{
    std::cout << "EvaluatorDoubleLep::Evaluate()\n";
    std::cout << m_chain->GetEntries() << "\n";
}