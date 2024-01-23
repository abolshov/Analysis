#ifndef M_TOOLS_HPP
#define M_TOOLS_HPP

#include <vector>

// finds closest daughters of particle at location part_idx in the event;
std::vector<int> GetNextGeneration(int part_idx, int const* mothers, int n_gen_part);

// finds all descendants of particle at location part_idx in the event;
// descendants are grouped into generations (inner vectors): generation 0 is particle itself, generation 1 are its immediate daughters, etc ...
std::vector<std::vector<int>> GetDescendants(int part_idx, int const* mothers, int n_gen_part);

// finds all particles by index that do not have daughters
std::vector<int> GetFinalParticles(int const* mothers, int n_gen_part);

#endif