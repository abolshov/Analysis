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
    if (tmp1.size() == 4)
    {
        for (auto const& idx: tmp1)
        {
            if (pdg_ids[idx] == B_ID) res[SIG::b] = idx;
            if (pdg_ids[idx] == BBAR_ID) res[SIG::bbar] = idx;
        }
    }

    int wplus = FindLast(WPLUS_ID, pdg_ids, n_gen_part);
    int wminus = FindLast(WMINUS_ID, pdg_ids, n_gen_part);
    auto w_plus_daughters = GetNextGeneration(wplus, mothers, n_gen_part);
    auto w_minus_daughters = GetNextGeneration(wminus, mothers, n_gen_part);

    if (w_plus_daughters.size() == 2)
    {
        if (std::abs(pdg_ids[w_plus_daughters[0]]) < B_ID && std::abs(pdg_ids[w_plus_daughters[1]]) < B_ID)
        {
            res[SIG::q1] = w_plus_daughters[0];
            res[SIG::q2] = w_plus_daughters[1];
        }
        else if (std::abs(pdg_ids[w_plus_daughters[0]]) > B_ID && std::abs(pdg_ids[w_plus_daughters[1]]) > B_ID)
        {
            res[SIG::l] = IsLepton(pdg_ids[w_plus_daughters[0]]) ? w_plus_daughters[0] : w_plus_daughters[1];
            res[SIG::nu] = IsNeutrino(pdg_ids[w_plus_daughters[1]]) ? w_plus_daughters[1] : w_plus_daughters[0];
        }
    }

    if (w_minus_daughters.size() == 2)
    {
        if (std::abs(pdg_ids[w_minus_daughters[0]]) < B_ID && std::abs(pdg_ids[w_minus_daughters[1]]) < B_ID)
        {
            res[SIG::q1] = w_minus_daughters[0];
            res[SIG::q2] = w_minus_daughters[1];
        }
        else if (std::abs(pdg_ids[w_minus_daughters[0]]) > B_ID && std::abs(pdg_ids[w_minus_daughters[1]]) > B_ID)
        {
            res[SIG::l] = IsLepton(pdg_ids[w_minus_daughters[0]]) ? w_minus_daughters[0] : w_minus_daughters[1];
            res[SIG::nu] = IsNeutrino(pdg_ids[w_minus_daughters[1]]) ? w_minus_daughters[1] : w_minus_daughters[0];
        }
    }

    if (std::any_of(res.begin(), res.end(), [](int x){ return x == -1; })) res.clear();
    return res;
}

bool CheckSignal(std::vector<int> const& signal, int const* mothers, int const* pdg_ids, int n_gen_part)
{
    if (signal.size() != N_SIG_PART)
    {
        return false;
    }
    
    int rad_idx = FindLast(RADION_ID, pdg_ids, n_gen_part);

    for (auto s: signal)
    {
        if (!IsDescOf(s, rad_idx, mothers))
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

int Match(int idx, KinematicData const kd_part, KinematicData kd_jet)
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

    if (best_match != hash.end()) return best_match->second;
    return -1;
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

std::unique_ptr<TH2F> EnergyMap(int const event_num, KinematicData const& kd, int const* mothers, int nbins, AxisRange xrange, AxisRange yrange)
{
    auto&& [xmin, xmax] = xrange;
    auto&& [ymin, ymax] = yrange;
    auto hist = std::make_unique<TH2F>(Form("evt_%d_EM", event_num), Form("Event %d", event_num), nbins, xmin, xmax, nbins, ymin, ymax); 
    int n_parts = kd.n;
    std::vector<int> final_parts = GetFinalParticles(mothers, n_parts);

    for (auto const& part_idx: final_parts)
    {
        TLorentzVector p4 = GetP4(kd, part_idx);
        hist->Fill(p4.Phi(), p4.Eta(), p4.E());
    }

    return hist;   
}

std::unique_ptr<TGraph> DRCone(KinematicData const& kd, int idx)
{
    TLorentzVector p4 = GetP4(kd, idx);
    double eta_0 = p4.Eta();
    double phi_0 = p4.Phi();

    std::vector<double> phis(N_POINTS + 1);
    std::vector<double> etas(N_POINTS + 1);

    for (int i = 0; i < N_POINTS + 1; ++i)
    {
        double t = 2*3.14/(N_POINTS)*i;
        phis[i] = phi_0 + DR_THRESH*std::cos(t);
        etas[i] = eta_0 + DR_THRESH*std::sin(t);
    }
    return std::make_unique<TGraph>(N_POINTS + 1, phis.data(), etas.data());
}

void DrawEventMap(MatchKinematics const& match_kin, MatchIndex const& match_index, int evt_num, std::pair<int*, int*> ptrs)
{
    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);

    auto const& [mothers, pdg_ids] = ptrs;
    auto const& [genpart, genjet] = match_kin;

    std::unique_ptr<TH2F> en_map = EnergyMap(evt_num, genpart, mothers);   
    en_map->SetStats(0);
    en_map->GetXaxis()->SetTitle("phi");
    en_map->GetYaxis()->SetTitle("eta");
    en_map->Draw("colz");

    std::vector<std::unique_ptr<TGraph>> part_cones, jet_cones;
    // cannot draw vectors in the loop as they will go out of scope after the loop is executed, thus need to save them
    for (auto const& match_idx_pair: match_index)
    {
        auto [part_idx, jet_idx] = match_idx_pair;
        auto&& part_cone = DRCone(genpart, part_idx);
        auto&& jet_cone = DRCone(genjet, jet_idx);

        part_cones.push_back(std::move(part_cone));
        jet_cones.push_back(std::move(jet_cone));
    }

    auto leg = std::make_unique<TLegend>(0.15, 0.1, 0.35, 0.3);
    for (int i = 0; i < static_cast<int>(part_cones.size()); ++i)
    {
        auto [part_idx, jet_idx] = match_index[i];
        part_cones[i]->SetLineWidth(2);
        part_cones[i]->SetLineColor(i+6);
        part_cones[i]->Draw("same");
        jet_cones[i]->SetLineWidth(2);
        jet_cones[i]->SetLineColor(i+6);
        jet_cones[i]->Draw("same");
        leg->AddEntry(part_cones[i].get(), Form("pdg ID = %d", pdg_ids[part_idx]));
    }
    leg->Draw();

    c1->SaveAs(Form("EnergyMaps/Event_%d.png", evt_num));
}

TLorentzVector Cone(int target, int const* mothers, KinematicData const& kd, double rad)
{
    int n_gen_part = kd.n;
    TLorentzVector target_p4 = GetP4(kd, target);
    std::vector<int> finals = GetFinalParticles(mothers, n_gen_part);
    TLorentzVector cone;
    for (auto idx: finals)
    {
        TLorentzVector p4 = GetP4(kd, idx);
        double dR = target_p4.DeltaR(p4);
        if (dR < rad)
        {
            bool belongs_to_decay = IsDescOf(idx, target, mothers);
            if (belongs_to_decay)
            {
                cone += p4;
            }
        }
    }
    return cone;
}

bool PassGenJetCut(std::vector<TLorentzVector> const& jets)
{
    for (auto const& j: jets)
    {
        if (j.Pt() < MIN_GENJET_PT || std::abs(j.Eta()) > MAX_GENJET_ETA)
        {
            return false;
        }
    }
    return true;
}

bool IsIsolatedLepton(TLorentzVector const& lep, std::vector<TLorentzVector> const& jets)
{
    for (auto const& j: jets)
    {
        if (j.DeltaR(lep) < DR_THRESH)
        {
            return false;
        }
    }
    return true;
}
