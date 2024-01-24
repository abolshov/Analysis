#ifndef M_TOOLS_HPP
#define M_TOOLS_HPP

#include <vector>
#include <iostream>

// finds closest daughters of particle at location part_idx in the event;
std::vector<int> GetNextGeneration(int part_idx, int const* mothers, int n_gen_part);

// finds all descendants of particle at location part_idx in the event;
// descendants are grouped into generations (inner vectors): generation 0 is particle itself, generation 1 are its immediate daughters, etc ...
std::vector<std::vector<int>> GetDescendants(int part_idx, int const* mothers, int n_gen_part);

// finds all particles by index that do not have daughters
std::vector<int> GetFinalParticles(int const* mothers, int n_gen_part);

// by default onlyh prints pdg_id(mother_idx) (in that format)
// if print_idx is true, prints pdg_id[idx](mother_idx)
void PrintDecay(std::vector<std::vector<int>> const& decay, int const* pdg_ids, int const* mothers, bool print_idx = false, std::ostream& stream = std::cout);

// returns index of the last occurence of particle with pdgId target_id
int FindLast(int target_id, int const* pdg_ids, int n_gen_part);
#endif