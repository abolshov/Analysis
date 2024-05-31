#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include <optional>
#include <map>
#include <vector>
#include <memory>

#include "EventData.hpp"
#include "Constants.hpp"
#include "SignalData.hpp"

class Selector 
{
    public:
    // fills maps 
    // head of the decay tree will be a particle with head_id
    void BuildTree(int head_id, std::unique_ptr<EventData> const& data);
    void PrintTree(std::unique_ptr<EventData> const& data) const;

    std::optional<SignalData> Select(std::unique_ptr<EventData> const& data) const;

    private:
    // computes indices of first-generation daughters of particle part_idx
    std::vector<int> GetNextGen(int part_idx, int const* mothers, int n_gen_part) const;
    std::vector<int> FindFirst(int target_id, int mother_idx, int const* pdg_ids, int const* mothers) const;
    std::vector<int> FindLast(int target_id, int mother_idx, int const* pdg_ids, int const* mothers) const;
    std::map<int, std::vector<int>> FindAllDesc(int target_id, int mother_idx, int const* pdg_ids, int const* mothers) const;
    // computes all first generation occurences and last generation occurences of particle target_id descending from mother_idx
    // if nothing found returns a pair of empty vectors
    // if first == last returns two copies of the same vector
    std::pair<std::vector<int>, std::vector<int>> FindFirstLast(int target_id, int mother_idx, int const* pdg_ids, int const* mothers) const;
    bool IsDescOf(int cand_idx, int parent_idx, int const* mothers) const;

    // map of generation to indices of particles in this generation
    // generation 0 will consist of the single particle - first particle head_id
    std::map<int, std::vector<int>> m_generations;
    // maps generation to list of indicating stable and unstable particle (stable = true, unstable = false)
    // stable = has no daughters 
    std::map<int, std::vector<bool>> m_stable;
};

#endif