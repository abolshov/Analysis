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

    if (std::any_of(res.begin(), res.end(), [](int x){ return x == -1; })) res.clear();
    return res;
}

bool CheckSignal(std::vector<int> const& signal, int const* mothers, int const* pdg_ids)
{
    if (signal.size() != N_SIG_PART) return false;
    std::vector<int> copy = signal;
    if (std::unique(copy.begin(), copy.end()) != copy.end()) return false;
    if (std::any_of(copy.begin(), copy.end(), [](int x){ return x == -1; })) return false;

    int h1 = signal[SIG::h1];
    int h2 = signal[SIG::h2];
    int b = signal[SIG::b];
    int bbar = signal[SIG::bbar];
    int q1 = signal[SIG::q1];
    int q2 = signal[SIG::q2];
    int l = signal[SIG::l];
    int nu = signal[SIG::nu];

    // std::cout << "h1 = " << h1 << "\n"
    //           << "h2 = " << h2 << "\n"
    //           << "b = " << b << "\n"
    //           << "bbar = " << bbar << "\n"
    //           << "q1 = " << q1 << "\n"
    //           << "q2 = " << q2 << "\n"
    //           << "l = " << l << "\n"
    //           << "nu = " << nu << "\n";
              

    bool b_quarks = (std::abs(pdg_ids[b]) == B_ID && std::abs(pdg_ids[bbar]) == B_ID);
    // std::cout << std::boolalpha;
    // std::cout << "b_quarks = " << b_quarks << "\n";
    if (!b_quarks) return false;
    
    bool h1tobb = (mothers[b] == h1 && mothers[bbar] == h1);
    bool h2tobb = (mothers[b] == h2 && mothers[bbar] == h2);

    // std::cout << "h1tobb = " << h1tobb << "\n" 
    //           << "h2tobb = " << h2tobb << "\n";

    if (!h1tobb && !h2tobb) return false;
    if (h1tobb && h2tobb) return false;

    int hadr_w = -1;
    if (mothers[q1] != mothers[q2])
    {
        return false;
    }
    else
    {
        hadr_w = mothers[q1];
    }

    // std::cout << "hadr_w = " << hadr_w << ", pid = " << pdg_ids[hadr_w] << "\n";

    if (std::abs(pdg_ids[hadr_w]) != W_ID) return false;

    int lep_w = -1;
    if (mothers[l] != mothers[nu])
    {
        return false;
    }
    else
    {
        lep_w = mothers[l];
    }

    // std::cout << "lep_w = " << lep_w << ", pid = " << pdg_ids[lep_w] << "\n";

    if (std::abs(pdg_ids[lep_w]) != W_ID) return false;

    if (mothers[lep_w] != mothers[hadr_w]) return false;
    int h_toWW = mothers[lep_w];

    // std::cout << "h_toWW = " << h_toWW << "\n";
    // std::cout << "pdg_ids[h_toWW] = " << pdg_ids[h_toWW] << "\n";
    if (pdg_ids[h_toWW] != HIGGS_ID) return false;

    return true;
}

TLorentzVector GetP4(KinematicData const& kd, int idx)
{
    TLorentzVector p;
    auto const& [pt, eta, phi, m, n] = kd;
    p.SetPtPhiEtaM(pt[idx], eta[idx], phi[idx], m[idx]);
    return p;
}