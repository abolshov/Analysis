#ifndef M_TOOLS_HPP
#define M_TOOLS_HPP

#include <vector>
#include <iostream>
#include <algorithm>
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

static constexpr int RADION_ID = 35;
static constexpr int HIGGS_ID = 25;
static constexpr int W_ID = 24;
static constexpr int TAU_ID = 15;
static constexpr int NU_TAU_ID = 16;
static constexpr int WPLUS_ID = 24;
static constexpr int WMINUS_ID = -24;
static constexpr int B_ID = 5;
static constexpr int BBAR_ID = -5;

static constexpr double DR_THRESH = 0.5;
static constexpr int N_POINTS = 20;

static constexpr double MIN_GENJET_PT = 20.0;
static constexpr double MAX_GENJET_ETA = 2.5;

static constexpr double MIN_LEP_PT = 5.0;
static constexpr double MAX_LEP_ETA = 2.5;

// specifies order of signal (hh->bbWW->bbqqlv) particles
enum SIG { X, H1, H2, b, bbar, LepWfirst, HadWfirst, HadWlast, q1, q2, LepWlast, l, nu };

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

// returns p4 of particle at index idx
TLorentzVector GetP4(KinematicData const& kd, int idx);

// matches gen jet to gen particle by dR
// returns index of jet matched to quark
// if matching fails returns -1
int Match(int idx, KinematicData const& kd_part, KinematicData const& kd_jet);

// computes min dR between all particles
double MinDeltaR(std::vector<TLorentzVector> const& parts);

// compute energy map of an event
using AxisRange = std::pair<double, double>;
std::unique_ptr<TH2F> EnergyMap(int const event_num, KinematicData const& kd, int const* mothers, int nbins = 30, AxisRange xrange = {-6.5, 6.5}, AxisRange yrange = {-6.5, 6.5});

// Make a graph with a marker centered at particle idx in (eta, phi) space and 0.4 dR circle around it
std::unique_ptr<TGraph> DRCone(KinematicData const& kd, int idx);

using MatchKinematics = std::pair<KinematicData, KinematicData>;
using MatchIndex = std::vector<std::pair<int, int>>; // contains pairs: {part_idx, jet_matched_to_part_idx}
void DrawEventMap(MatchKinematics const& match_kin, MatchIndex const& match_index, int evt_num, std::pair<int*, int*> ptrs);

// checks if all matched jets pass genjet acceptance cut
// returns true if all jets satisfy selection criteria
bool PassGenJetCut(std::vector<TLorentzVector> const& jets);

// checks if lepton passes acceptance cut
// returns true if lepton satisfies selection criteria
inline bool PassLeptonCut(TLorentzVector const& lep)
{
    return lep.Pt() > MIN_LEP_PT && std::abs(lep.Eta()) < MAX_LEP_ETA;
}

// checks that lepton is separated from jets by dR
bool IsIsolatedLepton(TLorentzVector const& lep, std::vector<TLorentzVector> const& jets);
bool IsIsolatedLepton(TLorentzVector const& lep, KinematicData const& kd);

// checks if matched objects preserve relationship between underlying objetcts (pt)
using MatchedPair = std::pair<TLorentzVector, TLorentzVector>;
inline bool ConsistentMatch(MatchedPair const& mp1, MatchedPair const& mp2)
{
    auto const& [q1, j1] = mp1;
    auto const& [q2, j2] = mp2;

    return (q1.Pt() > q2.Pt() ? j1.Pt() > j2.Pt() : j1.Pt() < j2.Pt());
};

// returns indices of jets in acceptance region (pt > 25 and |eta| < 2.5)
std::vector<int> GetAcceptJets(KinematicData const& kd);

#endif