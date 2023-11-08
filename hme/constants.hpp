#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "tools.hpp"

constexpr float MET_SIGMA = 25.2f;
constexpr size_t NUM_ITER = 1000;
constexpr RescFact FAIL = {-1.0f, -1.0f};
// constexpr float ONSHELL_WEIGHT = 0.15f; // 0.245 good for 350 GEV
constexpr float ONSHELL_WEIGHT = 0.5f;
constexpr float LOWER_BOUND = 50.0f;
constexpr float UPPER_BOUND = 70.0f;
constexpr float MASS_LIMIT = 90.0f;
constexpr size_t NUM_RESC_INTEGRATIONS = 1000;
constexpr size_t NUM_MASS_INTEGRATIONS = 100;
constexpr float W_MASS_THRESHOLD = 62.0f;

#endif