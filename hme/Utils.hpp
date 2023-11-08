#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <functional>

#include "TLorentzVector.h"

extern TLorentzVector const zero;

enum PART {b1, b2, q1, q2, l, met};

bool ValidEvent(std::vector<TLorentzVector> const& particles);
bool PartonCut(std::vector<TLorentzVector> const& particles); 
bool GenJetCut(std::vector<TLorentzVector> const& particles);
void Print(TLorentzVector const& p, bool EXYZ);

#endif