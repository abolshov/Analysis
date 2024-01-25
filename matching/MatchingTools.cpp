#include "MatchingTools.hpp"

std::vector<int> GetNextGeneration(int part_idx, int const* mothers, int n_gen_part)
{
    std::vector<int> gen;
    for (int i = 0; i < n_gen_part; ++i)
    {
        if (mothers[i] == part_idx) gen.push_back(i);
    }
    if (part_idx == -1) gen.clear();
    return gen;
}

std::vector<std::vector<int>> GetDescendants(int part_idx, int const* mothers, int n_gen_part)
{
    std::vector<std::vector<int>> descendants;
    descendants.push_back({part_idx});
    // int gener_counter = 1;
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
        // ++gener_counter;
        cur_gen = next_gen;
    }
    if (part_idx == -1) descendants.clear();
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

void PrintDecay(std::vector<std::vector<int>> const& decay, int const* pdg_ids, int const* mothers, bool print_idx, std::ostream& stream)
{
    for (int gen = 0; gen < static_cast<int>(decay.size()); ++gen)
    {
        stream << "generation " << gen << ":\n";
        stream << "\t";
        auto cur_gen = decay[gen];
        for (auto idx: cur_gen)
        {
            if (print_idx)
            {
                stream << pdg_ids[idx] << "[" << idx << "]" << "(" << mothers[idx] << ") ";
            }
            else
            {
                stream << pdg_ids[idx] << "(" << mothers[idx] << ") ";
            }
            
        }
        stream << "\n";
    }
}

int FindLast(int target_id, int const* pdg_ids, int n_gen_part)
{
    int res = -1;
    for (int i = n_gen_part - 1; i > -1; --i)
    {
        if (pdg_ids[i] == target_id)
        {
            res = i;
            break;
        }
    }
    return res;
}

std::vector<int> FindSpecificDescendants(std::vector<int> const& desc_range, int mother_idx, int const* mothers, int const* pdg_ids, int n_gen_part)
{
    std::vector<int> res;
    auto descendants = GetNextGeneration(mother_idx, mothers, n_gen_part);
    for (auto const& desc_pdg_id: desc_range)
    {
        for (auto const& idx: descendants)
        {
            if (desc_pdg_id == pdg_ids[idx]) res.push_back(idx);
        }
    }
    return res;
}

std::vector<int> GetSignal(int const* pdg_ids, int const* mothers, int n_gen_part)
{
    std::vector<int> res(N_SIG_PART, -1);
    // res.reserve(N_SIG_PART);
    int rad_idx = FindLast(RADION_ID, pdg_ids, n_gen_part);
    
    // auto gen = GetNextGeneration(rad_idx, mothers, n_gen_part);
    // tmp vector to store generation being examined currently
    auto tmp1 = FindSpecificDescendants({HIGGS_ID}, rad_idx, mothers, pdg_ids, n_gen_part);
    // H0->hh only
    if (tmp1.size() == 2)
    {
        res[SIG::h1] = tmp1[0];
        res[SIG::h2] = tmp1[1];
    }

    tmp1 = GetNextGeneration(res[SIG::h1], mothers, n_gen_part); 
    auto tmp2 = GetNextGeneration(res[SIG::h2], mothers, n_gen_part); 
    tmp1.insert(tmp1.end(), tmp2.begin(), tmp2.end()); // now tmp1 contains indices of all daughters of both higgses or is empty;
    
    // hh->bbWW only!
    int wplus = -1;
    int wminus = -1;
    if (tmp1.size() == 4)
    {
        for (auto const& idx: tmp1)
        {
            if (pdg_ids[idx] == B_ID) res[SIG::b] = idx;
            if (pdg_ids[idx] == BBAR_ID) res[SIG::bbar] = idx;
            if (pdg_ids[idx] == WPLUS_ID) wplus = idx;
            if (pdg_ids[idx] == WMINUS_ID) wminus = idx; 
        }
    }

    tmp1 = GetNextGeneration(wplus, mothers, n_gen_part); 
    tmp2 = GetNextGeneration(wminus, mothers, n_gen_part); 
    tmp1.insert(tmp1.end(), tmp2.begin(), tmp2.end()); // now tmp1 contains all daughters of both W

    if (tmp1.size() == 4)
    {
        // after sorting by abs(pdgIds) indices of quarks will always be the first two elements 
        // third elem will laways be lepton as std::abs(lep(pdgId)) < std::abs(nu(pdgId))
        std::sort(tmp1.begin(), tmp1.end(), [&pdg_ids](int x, int y){ return std::abs(pdg_ids[x]) < std::abs(pdg_ids[y]); });

        res[SIG::q1] = IsLightQuark(pdg_ids[tmp1[0]]) ? tmp1[0] : -1;
        res[SIG::q2] = IsLightQuark(pdg_ids[tmp1[1]]) ? tmp1[1] : -1;
        res[SIG::l] = IsLepton(pdg_ids[tmp1[2]]) ? tmp1[2] : -1;
        res[SIG::nu] = IsNeutrino(pdg_ids[tmp1[3]]) ? tmp1[3] : -1;
    }

    std::cout << "before the last check:\n";
    for (auto const& s: res)
    {
        std::cout << s << " "; 
    }
    std::cout << "\n";

    if (std::any_of(res.begin(), res.end(), [](int x){ return x == -1; })) res.clear();
    return res;
}