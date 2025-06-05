#include "EvaluatorDoubleLep.hpp"

#include <iostream>

void EvaluatorDoubleLep::Evaluate()
{
    ULong64_t n_events = m_chain->GetEntries();
    for (ULong64_t i = 0; i < n_events; ++i)
    {
        m_chain->GetEntry(i);
        // do stuff here
    }
}