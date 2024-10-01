#ifndef M_TOOLS_HPP
#define M_TOOLS_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>

#include "TLorentzVector.h"
#include "TH2.h"
#include "TGraph.h"

// signal does not include tau leptons and neutrinos
static const std::vector<int> SIG_LIGHT_QUARKS = {1, 2, 3, 4};
static const std::vector<int> SIG_LEPTONS = {11, 13};
static const std::vector<int> SIG_NEUTRINOS = {12, 14};

static const std::vector<int> LEPTONS = {11, 13, 15};
static const std::vector<int> NEUTRINOS = {12, 14, 16};

static constexpr int N_SIG_PART = 13;
static constexpr int N_JETS_AWAY_FROM_LEP = 4;

static constexpr int RADION_ID = 35;
static constexpr int HIGGS_ID = 25;
static constexpr int W_ID = 24;
static constexpr int TAU_ID = 15;
static constexpr int NU_TAU_ID = 16;
static constexpr int WPLUS_ID = 24;
static constexpr int WMINUS_ID = -24;
static constexpr int B_ID = 5;
static constexpr int BBAR_ID = -5;

static constexpr double DR_THRESH = 0.4;
static constexpr int N_POINTS = 20;

static constexpr double MIN_GENJET_PT = 20.0;
static constexpr double MAX_GENJET_ETA = 2.5;
static constexpr double PRIMARY_GENJET_PT = 20.0;
static constexpr double PRIMARY_GENJET_ETA = 2.5;

static constexpr double MIN_LEP_PT = 5.0;
static constexpr double MAX_LEP_ETA = 2.5;

// specifies order of signal (hh->bbWW->bbqqlv) particles
enum SIG { X, H_bb, H_WW, b, bbar, LepWfirst, HadWfirst, HadWlast, q1, q2, LepWlast, l, nu };

void Print(TLorentzVector const& p, bool EXYZ = false);

// finds closest daughters of particle at location part_idx in the event;
// returns indices of what found
std::vector<int> GetNextGeneration(int part_idx, int const* mothers, int n_gen_part);

// finds all descendants of particle at location part_idx in the event;
// descendants are grouped into generations (inner vectors): generation 0 is particle itself, generation 1 are its immediate daughters, etc ...
std::vector<std::vector<int>> GetDescendants(int part_idx, int const* mothers, int n_gen_part);

// finds all particles by index that do not have daughters
std::vector<int> GetFinalParticles(int const* mothers, int n_gen_part);

// finds stable descendants (that don't have any daughters) of particle at part_idx
std::vector<int> GetStableDescendants(int part_idx, int const* mothers, int n_gen_part);

// by default onlyh prints pdg_id(mother_idx) (in that format)
// if print_idx is true, prints pdg_id[idx](mother_idx)
void PrintDecay(std::vector<std::vector<int>> const& decay, int const* pdg_ids, int const* mothers, bool print_idx = false, std::ostream& stream = std::cout);

// returns index of the last occurence of particle with pdgId target_id
int FindLast(int target_id, int const* pdg_ids, int n_gen_part);

// returns indices of signal (see enum SIG) particles
// access specific particles via SIG::_
// head is position of the last heavy higgs, i.e. index where the decay starts 
// if any signal particle was not found returns empty vector
std::vector<int> GetSignal(int const* pdg_ids, int const* mothers, int n_gen_part);

// finds daughters with pdg ids in range desc_range of mother at mother_idx 
std::vector<int> FindSpecificDescendants(std::vector<int> const& desc_range, int mother_idx, int const* mothers, int const* pdg_ids, int n_gen_part);

// self-explanatory helper functions
inline bool IsSigLightQuark(int pdg_id) { return std::find(SIG_LIGHT_QUARKS.begin(), SIG_LIGHT_QUARKS.end(), std::abs(pdg_id)) != SIG_LIGHT_QUARKS.end(); }
inline bool IsSigLep(int pdg_id) { return std::find(SIG_LEPTONS.begin(), SIG_LEPTONS.end(), std::abs(pdg_id)) != SIG_LEPTONS.end(); }
inline bool IsSigNu(int pdg_id) { return std::find(SIG_NEUTRINOS.begin(), SIG_NEUTRINOS.end(), std::abs(pdg_id)) != SIG_NEUTRINOS.end(); }

inline bool IsNu(int pdg_id) { return std::find(NEUTRINOS.begin(), NEUTRINOS.end(), std::abs(pdg_id)) != NEUTRINOS.end(); }
inline bool IsLep(int pdg_id) { return std::find(LEPTONS.begin(), LEPTONS.end(), std::abs(pdg_id)) != LEPTONS.end(); }

// check validity of signal particles
bool HasOnlyEleMu(std::vector<int> const& signal, int const* mothers, int const* pdg_ids);

// checks of particle cand_idx is daughter of particle parent_idx
bool IsDescOf(int cand_idx, int parent_idx, int const* mothers);

struct KinematicData
{
    Float_t* pt = nullptr;
    Float_t* eta = nullptr;
    Float_t* phi = nullptr;
    Float_t* m = nullptr;
    int n = 0;
};

struct GenJetFlavorInfo
{
    Int_t* parton_flavor = nullptr;
    UChar_t* hadron_flavor = nullptr;
};

struct GenPartIdInfo
{
    Int_t* pdg_id = nullptr;
    Int_t* mothers = nullptr;
};

struct GenJetData
{
    KinematicData kd;
    GenJetFlavorInfo fi;
};

// returns p4 of particle at index idx
TLorentzVector GetP4(KinematicData const& kd, int idx);

// matches gen jet to gen particle by dR
// returns index of jet matched to quark
// if matching fails returns -1
int Match(int idx, KinematicData const& kd_part, KinematicData const& kd_jet);

// returns index in array jets
// no match == returns -1
int Match(TLorentzVector const& quark, std::vector<TLorentzVector> const& jets);

// if one match exists uses it to find best match by scanning all remaining jets and pick the one s.t. mass of the pair is closest to mass of the target particle
int FindSecondMatch(int first_match, TLorentzVector const& target, std::vector<TLorentzVector> const& jets);

// returns all posible matches (within 0.4 dR cone)
std::vector<int> Matches(int idx, KinematicData const& kd_part, KinematicData const& kd_jet);

TLorentzVector MatchDijet(TLorentzVector const& hadW, std::vector<TLorentzVector> const& jets);

// computes min dR between all particles
double MinDeltaR(std::vector<TLorentzVector> const& parts);
double MinDeltaR(TLorentzVector const& v, std::vector<TLorentzVector> const& other);

// checks that among passed jets there are at least 4 jets such that the lepton is separated from them by at least dR = 0.4
bool IsIsolatedLepton(TLorentzVector const& lep, std::vector<TLorentzVector> const& jets);
bool IsIsolatedLepton(TLorentzVector const& lep, KinematicData const& kd);

// performs selection of jets with pt > 10 and |eta| < 5
std::vector<int> PrimaryJetSelection(KinematicData const& kd);

inline bool PassLeptonCut(TLorentzVector const& lep)
{
    return lep.Pt() > MIN_LEP_PT && std::abs(lep.Eta()) < MAX_LEP_ETA;
}

template <typename Func>
std::pair<int, int> ChooseBestPair(std::vector<TLorentzVector> const& jets, Func func)
{
    size_t sz = jets.size();
    if (sz == 2)
    {
        return {0, 1};
    }


    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(sz*(sz-1)/2);
    for (size_t i = 0; i < sz; ++i)
    {
        for (size_t j = 0; j < sz; ++j)
        {
            pairs.emplace_back(i, j);
        }
    }

    double min_metric = 10e9;
    std::pair<int, int> res{-1, -1};
    for (auto const& p: pairs)
    {
        auto [i1, i2] = p;
        TLorentzVector const& jet1 = jets[i1];
        TLorentzVector const& jet2 = jets[i2];

        double metric = func(jet1, jet2);
        if (metric < min_metric)
        {
            min_metric = metric;
            res = p;
        }
    }

    return res;
}

template <typename Func>
std::pair<double, double> CalcJetPairStats(std::vector<TLorentzVector> const& jets, Func func)
{
    size_t n_jets = jets.size();
    if (n_jets == 2)
    {
        return {func(jets[0], jets[1]), 0.0};
    }

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(n_jets*(n_jets-1)/2);
    for (size_t i = 0; i < n_jets; ++i)
    {
        for (size_t j = 0; j < n_jets; ++j)
        {
            pairs.emplace_back(i, j);
        }
    }

    size_t sz = pairs.size();

    auto Add = [&jets, &func](double acc, std::pair<int, int> const& p)
    {
        auto [i1, i2] = p; 
        return acc + func(jets[i1], jets[i2]); 
    };
    double mean = std::accumulate(pairs.begin(), pairs.end(), 0.0, Add)/sz;

    auto Variance = [&mean, &sz, &jets, &func](double var, std::pair<int, int> const& p)
    {
        auto [i1, i2] = p; 
        double value = func(jets[i1], jets[i2]); 
        return var + ((value - mean)*(value - mean)/(sz - 1)); 
    };
    double sdtdev = std::sqrt(std::accumulate(pairs.begin(), pairs.end(), 0.0, Variance));
    return {mean, sdtdev};
}

#endif