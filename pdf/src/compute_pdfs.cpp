#include "TROOT.h"

#include <fstream>

#include "Definitions.hpp"
#include "Constants.hpp"
#include "EvaluatorDoubleLep.hpp"
#include "EvaluatorBase.hpp"

int main()
{
    TH1::AddDirectory(false);
    TH2::AddDirectory(false);

    std::ifstream files_dl("files_dl.txt");
    std::unique_ptr<EvaluatorBase> evaluator{std::make_unique<EvaluatorDoubleLep>(files_dl, Channel::DL)};
    evaluator->Evaluate();

    return 0;
}