#include "MatchingTools.hpp"

std::vector<int> GetNextGeneration(int part_idx, int const* mothers, int n_gen_part)
{
    std::vector<int> gen;
    for (int i = 0; i < n_gen_part; ++i)
    {
        if (mothers[i] == part_idx) gen.push_back(i);
    }
    return gen;
}

std::vector<std::vector<int>> GetDescendants(int part_idx, int const* mothers, int n_gen_part)
{
    std::vector<std::vector<int>> descendants;
    descendants.push_back({part_idx});
    int gener_counter = 1;
    auto cur_gen = GetNextGeneration(part_idx, mothers, n_gen_part);
    while (!cur_gen.empty())
    {
        descendants.push_back(cur_gen);
        std::vector<int> next_gen;
        for (auto const& part: cur_gen)
        {
            auto part_daughters = GetNextGeneration(part, mothers, n_gen_part);
            next_gen.insert(next_gen.end(), part_daughters.begin(), part_daughters.end());
        }
        ++gener_counter;
        cur_gen = next_gen;
    }
    return descendants;
}

std::vector<int> GetFinalParticles(int const* mothers, int n_gen_part)
{
    std::vector<int> finals;
    for (int i = 0; i < n_gen_part; ++i)
    {
        if (GetNextGeneration(i, mothers, n_gen_part).empty()) finals.push_back(i);
    }
    return finals;
}