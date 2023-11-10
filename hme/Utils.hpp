#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <functional>

#include "TLorentzVector.h"

TLorentzVector const zero(0.0f, 0.0f, 0.0f, 0.0f);

enum PART {b1, b2, q1, q2, l, met};

bool ValidEvent(std::vector<TLorentzVector> const& particles);
bool PartonCut(std::vector<TLorentzVector> const& particles); 
bool GenJetCut(std::vector<TLorentzVector> const& particles);
void Print(TLorentzVector const& p, bool EXYZ);

bool HasZeroParticle(std::vector<TLorentzVector> const& particles);
bool IsIdenticalPair(TLorentzVector const& p1, TLorentzVector const& p2);

#endif