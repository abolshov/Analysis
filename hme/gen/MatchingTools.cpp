#include "MatchingTools.hpp"

#include <map>

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"

void Print(TLorentzVector const& p, bool EXYZ)
{
    if (EXYZ) 
    {
        std::cout << "(" << p.E() << ", " << p.X() << ", " << p.Y() << ", " << p.Z() << ")\n";
    }
    else
    {
        std::cout << "(" << p.Pt() << ", " << p.Eta() << ", " << p.Phi() << ", " << p.M() << ")\n";
    }
}

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
        if (GetNextGeneration(i, mothers, n_gen_part).empty()) 
            finals.push_back(i);
    }
    return finals;
}

std::vector<int> GetStableDescendants(int part_idx, int const* mothers, int n_gen_part)
{
    std::vector<int> res;
    auto stable = GetFinalParticles(mothers, n_gen_part);
    for (auto idx: stable)
    {
        if (IsDescOf(idx, part_idx, mothers))
        {
            res.push_back(idx);
        }
    }
    return res;
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
    int rad_idx = FindLast(RADION_ID, pdg_ids, n_gen_part);
    res[SIG::X] = rad_idx;
    
    // tmp vector to store generation being examined currently
    auto tmp1 = FindSpecificDescendants({HIGGS_ID}, rad_idx, mothers, pdg_ids, n_gen_part);
    // H0->hh only
    if (tmp1.size() == 2)
    {
        res[SIG::H_bb] = tmp1[0];
        res[SIG::H_WW] = tmp1[1];
    }

    tmp1 = GetNextGeneration(res[SIG::H_bb], mothers, n_gen_part); 
    auto tmp2 = GetNextGeneration(res[SIG::H_WW], mothers, n_gen_part); 

    tmp1.insert(tmp1.end(), tmp2.begin(), tmp2.end()); // now tmp1 contains indices of all daughters of both higgses or is empty;
    
    // hh->bbWW only!
    int first_w_plus = -1;
    int first_w_minus = -1;
    if (tmp1.size() == 4)
    {
        if (std::abs(pdg_ids[tmp1[0]]) != B_ID)
        {
            std::swap(res[SIG::H_bb], res[SIG::H_WW]);
        }

        for (auto const& idx: tmp1)
        {
            if (pdg_ids[idx] == B_ID)
            {
                res[SIG::b] = idx;
            }
            else if (pdg_ids[idx] == BBAR_ID) 
            {
                res[SIG::bbar] = idx;
            }
            else if (pdg_ids[idx] == WPLUS_ID)
            {
                first_w_plus = idx;
            }
            else if (pdg_ids[idx] == WMINUS_ID)
            {
                first_w_minus = idx;
            }
        }
    }

    int wplus = FindLast(WPLUS_ID, pdg_ids, n_gen_part);
    int wminus = FindLast(WMINUS_ID, pdg_ids, n_gen_part);
    auto w_plus_daughters = GetNextGeneration(wplus, mothers, n_gen_part);
    auto w_minus_daughters = GetNextGeneration(wminus, mothers, n_gen_part);

    if (w_plus_daughters.size() == 2)
    {
        if (std::abs(pdg_ids[w_plus_daughters[0]]) <= B_ID && std::abs(pdg_ids[w_plus_daughters[1]]) <= B_ID)
        {
            res[SIG::q1] = w_plus_daughters[0];
            res[SIG::q2] = w_plus_daughters[1];
            res[SIG::HadWlast] = wplus;
            res[SIG::HadWfirst] = (IsDescOf(wplus, first_w_plus, mothers) || wplus == first_w_plus) ? first_w_plus : first_w_minus;
        }
        else if (std::abs(pdg_ids[w_plus_daughters[0]]) > B_ID && std::abs(pdg_ids[w_plus_daughters[1]]) > B_ID)
        {
            // res[SIG::l] = IsSigLep(pdg_ids[w_plus_daughters[0]]) ? w_plus_daughters[0] : w_plus_daughters[1];
            // res[SIG::nu] = IsSigNu(pdg_ids[w_plus_daughters[1]]) ? w_plus_daughters[1] : w_plus_daughters[0];
            res[SIG::l] = IsLep(pdg_ids[w_plus_daughters[0]]) ? w_plus_daughters[0] : w_plus_daughters[1];
            res[SIG::nu] = IsNu(pdg_ids[w_plus_daughters[1]]) ? w_plus_daughters[1] : w_plus_daughters[0];
            res[SIG::LepWlast] = wplus;
            res[SIG::HadWfirst] = first_w_plus;
            res[SIG::LepWfirst] = (IsDescOf(wplus, first_w_plus, mothers) || wplus == first_w_plus) ? first_w_plus : first_w_minus;
        }
    }

    if (w_minus_daughters.size() == 2)
    {
        if (std::abs(pdg_ids[w_minus_daughters[0]]) <= B_ID && std::abs(pdg_ids[w_minus_daughters[1]]) <= B_ID)
        {
            res[SIG::q1] = w_minus_daughters[0];
            res[SIG::q2] = w_minus_daughters[1];
            res[SIG::HadWlast] = wminus;
            res[SIG::HadWfirst] = (IsDescOf(wminus, first_w_minus, mothers) || wminus == first_w_minus) ? first_w_minus : first_w_plus;
        }
        else if (std::abs(pdg_ids[w_minus_daughters[0]]) > B_ID && std::abs(pdg_ids[w_minus_daughters[1]]) > B_ID)
        {
            // res[SIG::l] = IsSigLep(pdg_ids[w_minus_daughters[0]]) ? w_minus_daughters[0] : w_minus_daughters[1];
            // res[SIG::nu] = IsSigNu(pdg_ids[w_minus_daughters[1]]) ? w_minus_daughters[1] : w_minus_daughters[0];
            res[SIG::l] = IsLep(pdg_ids[w_minus_daughters[0]]) ? w_minus_daughters[0] : w_minus_daughters[1];
            res[SIG::nu] = IsNu(pdg_ids[w_minus_daughters[1]]) ? w_minus_daughters[1] : w_minus_daughters[0];
            res[SIG::LepWlast] = wminus;
            res[SIG::LepWfirst] = (IsDescOf(wminus, first_w_minus, mothers) || wminus == first_w_minus) ? first_w_minus : first_w_plus;
        }
    }

    if (std::any_of(res.begin(), res.end(), [](int x){ return x == -1; })) 
    {
        res.clear();
    }
    return res;
}

bool HasOnlyEleMu(std::vector<int> const& signal, int const* mothers, int const* pdg_ids)
{
    if (signal.size() != N_SIG_PART)
    {
        return false;
    }
    
    for (auto s: signal)
    {
        if (pdg_ids[s] == RADION_ID)
        {
            continue;
        }
        if (!IsDescOf(s, signal[SIG::X], mothers))
        {
            return false;
        }
    }

    if (std::abs(signal[SIG::l]) == TAU_ID || std::abs(signal[SIG::nu]) == NU_TAU_ID)
    {
        return false;
    }

    return true;
}

bool IsDescOf(int cand_idx, int parent_idx, int const* mothers)
{
    int mother_idx = mothers[cand_idx];
    while(mother_idx != -1)
    {
        if (mother_idx == parent_idx) return true;
        int new_mother_idx = mothers[mother_idx];
        mother_idx = new_mother_idx;
    }
    return false;
}

TLorentzVector GetP4(KinematicData const& kd, int idx)
{
    TLorentzVector p;
    auto const& [pt, eta, phi, m, n] = kd;
    p.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], m[idx]);
    return p;
}

int Match(int idx, KinematicData const& kd_part, KinematicData const& kd_jet)
{
    TLorentzVector ppart = GetP4(kd_part, idx);
    std::map<double, int> hash;
    int n_jets = kd_jet.n;
    for (int i = 0; i < n_jets; ++i)
    {
        TLorentzVector pjet = GetP4(kd_jet, i);
        double dr = ppart.DeltaR(pjet);
        if (dr < DR_THRESH)
        {
            hash.insert({dr, i});
        }
    }

    auto best_match = hash.begin();

    if (best_match != hash.end())
    {
        return best_match->second;
    }
    return -1;
}

int Match(TLorentzVector const& quark, std::vector<TLorentzVector> const& jets)
{
    std::map<double, int> hash;
    int n_jets = jets.size();
    for (int i = 0; i < n_jets; ++i)
    {
        TLorentzVector const& jet = jets[i];
        double dr = quark.DeltaR(jet);
        if (dr < DR_THRESH)
        {
            hash.insert({dr, i});
        }
    }

    auto best_match = hash.begin();

    if (best_match != hash.end())
    {
        return best_match->second;
    }
    return -1;
}

int FindSecondMatch(int first_match, TLorentzVector const& target, std::vector<TLorentzVector> const& jets)
{
    TLorentzVector const& jet1 = jets[first_match];
    int res = -1;
    double min_dm = 10e9;
    int sz = jets.size();
    for (int i = 0; i < sz; ++i)
    {
        if (i == first_match)
        {
            continue;
        }
        TLorentzVector const& jet2 = jets[i];
        double dm = std::abs(target.M() - (jet1 + jet2).M());
        if (dm < min_dm)
        {
            min_dm = dm;
            res = i;
        }
    }
    return res; 
}

std::vector<int> Matches(int idx, KinematicData const& kd_part, KinematicData const& kd_jet)
{
    std::vector<int> matches;
    TLorentzVector p = GetP4(kd_part, idx);
    int n_jets = kd_jet.n;
    for (int i = 0; i < n_jets; ++i)
    {
        TLorentzVector j = GetP4(kd_jet, i);
        double dr = p.DeltaR(j);
        if (dr < DR_THRESH)
        {
            matches.push_back(i);
        }
    }
    return matches;
}

double MinDeltaR(std::vector<TLorentzVector> const& parts)
{
    int n = parts.size();
    double min_dR = parts[0].DeltaR(parts[1]);
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double dr = parts[i].DeltaR(parts[j]);
            if (dr < min_dR)
            {
                min_dR = dr;
            }
        }
    }
    return min_dR;
}

TLorentzVector MatchDijet(TLorentzVector const& hadW, std::vector<TLorentzVector> const& jets)
{
    size_t sz = jets.size();
    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(sz*(sz-1)/2);
    for (size_t i = 0; i < sz; ++i)
    {
        for (size_t j = 0; j < sz; ++j)
        {
            pairs.emplace_back(i, j);
        }
    }

    // double min_dm = 10e9;
    // double min_dr = 6.0;
    // TLorentzVector res{};
    // for (auto const& p: pairs)
    // {
    //     auto [i1, i2] = p;
    //     TLorentzVector dijet = jets[i1] + jets[i2];
    //     double dr = hadW.DeltaR(dijet);

    //     if (dr < DR_THRESH && dr < min_dr)
    //     {
    //         min_dr = dr;
    //         double dm = std::abs(dijet.M() - hadW.M());
    //         if (dm < min_dm)
    //         {
    //             res = dijet;
    //         }
    //     }
    // }

    // return res;

    std::vector<TLorentzVector> close_dijets;
    for (auto const& p: pairs)
    {
        auto [i1, i2] = p;
        TLorentzVector dijet = jets[i1] + jets[i2];
        double dr = hadW.DeltaR(dijet);

        if (dr < DR_THRESH)
        {
            close_dijets.push_back(dijet);
        }
    }

    // if (close_dijets.empty())
    // {
    //     return TLorentzVector{};
    // }

    // auto Cmp = [&hadW](TLorentzVector const& p1, TLorentzVector const& p2)
    // {
    //     return hadW.DeltaR(p1) < hadW.DeltaR(p2);
    // };
    // std::sort(close_dijets.begin(), close_dijets.end(), Cmp);
    // return *close_dijets.begin();

    auto Comparator = [&hadW](TLorentzVector const& p1, TLorentzVector const& p2)
    { 
        return std::abs(p1.M() - hadW.M()) < std::abs(p2.M() - hadW.M()); 
    };

    if (close_dijets.empty())
    {
        auto copy = jets;
        std::sort(copy.begin(), copy.end(), Comparator);
        return *copy.begin();
    }

    std::sort(close_dijets.begin(), close_dijets.end(), Comparator);
    return *close_dijets.begin();
}


bool IsIsolatedLepton(TLorentzVector const& lep, std::vector<TLorentzVector> const& jets)
{
    int count = jets.size();
    for (auto const& j: jets)
    {
        if (j.DeltaR(lep) < DR_THRESH)
        {
            --count;
        }
    }
    return count >= N_JETS_AWAY_FROM_LEP;
}

bool IsIsolatedLepton(TLorentzVector const& lep, KinematicData const& kd)
{
    auto const& [pt, eta, phi, m, n] = kd;
    int count = n;
    for (int j = 0; j < n; ++j)
    {
        TLorentzVector p = GetP4(kd, j);
        if (p.DeltaR(lep) < DR_THRESH)
        {
            --count;
        }
    }
    return count >= N_JETS_AWAY_FROM_LEP;
}

std::vector<int> PrimaryJetSelection(KinematicData const& kd)
{
    auto const& [pt, eta, phi, m, n] = kd;
    std::vector<int> res;
    for (int i = 0; i < n; ++i)
    {
        if (pt[i] > PRIMARY_GENJET_PT && std::abs(eta[i]) < PRIMARY_GENJET_ETA)
        {
            res.push_back(i);
        }
    }
    return res;
}

double MinDeltaR(TLorentzVector const& v, std::vector<TLorentzVector> const& other)
{
    double res = 10.0;
    for (auto const& p: other)
    {
        double dr = v.DeltaR(p);
        res = std::min(dr, res);
    }
    return res;
}
