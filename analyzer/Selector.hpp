#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include <optional>
#include <map>

// #include "SignalData.hpp"

class Selector 
{
    public:
    // Selector(std::unique_ptr<EventData> const& data);

    // fills maps
    void BuildTree(std::unique_ptr<EventData> const& data);
    void PrintTree(std::unique_ptr<EventData> const& data) const;

    private:
    // computes indices of first-generation daughters of particle part_idx
    std::vector<int> GetNextGen(int part_idx, int const* mothers, int n_gen_part) const;

    // map of generation to indices of particles in this generation
    // generation 0 will consist of the single particle - first produced Radion (X)
    std::map<int, std::vector<int>> m_generations;
    // maps generation to list of indicating stable and unstable particle (stable = true, unstable = false)
    // stable = has no daughters 
    std::map<int, std::vector<bool>> m_stable;
    // std::unique_ptr<EventData> const& m_data;
};

#endif