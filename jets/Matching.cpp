#include <iostream>
#include <algorithm>
#include <set>
#include <vector>

#include "TLorentzVector.h"

#include "Matching.hpp"
#include "Utils.hpp"

GenPartIndex IsDiHiggsSL(Int_t const* pdgIds, Int_t const* motherIdxs, Int_t const* statuses, Int_t nGenPart)
{
    GenPartIndex res;
    auto& [h1, h2, w1, w2, b1, b2, q1, q2, l, nu] = res;
    Int_t const* pdgIds_begin = pdgIds;
    Int_t const* pdgIds_end = pdgIds + nGenPart;

    auto p1 = std::find(pdgIds_begin, pdgIds_end, 25);
    auto p2 = std::find(p1 + 1, pdgIds_end, 25);

    Int_t idx1, idx2;
    idx1 = p1 - pdgIds_begin;
    idx2 = p2 - pdgIds_begin;

    if (motherIdxs[idx1] == motherIdxs[idx2])
    {
        h1 = idx1;
        h2 = idx2;
    }

    auto IsW = [](Int_t pdgId) { return std::abs(pdgId) == 24; };
    p1 = std::find_if(pdgIds_begin, pdgIds_end, IsW);
    p2 = std::find_if(p1 + 1, pdgIds_end, IsW);

    idx1 = p1 - pdgIds_begin;
    idx2 = p2 - pdgIds_begin;

    if (motherIdxs[idx1] == motherIdxs[idx2])
    {
        w1 = idx1;
        w2 = idx2;
    }
    
    auto IsBquark = [](Int_t pdgId) { return std::abs(pdgId) == 5; };
    p1 = std::find_if(pdgIds_begin, pdgIds_end, IsBquark);
    p2 = std::find_if(p1 + 1, pdgIds_end, IsBquark);

    idx1 = p1 - pdgIds_begin;
    idx2 = p2 - pdgIds_begin;

    if (motherIdxs[idx1] == motherIdxs[idx2] && statuses[idx1] == 23 && statuses[idx2] == 23)
    {
        b1 = idx1;
        b2 = idx2;
    }

    std::set<Int_t> light_quarks = {1, 2, 3, 4};
    auto IsLightQuark = [&light_quarks](Int_t pdgId)
    {
        auto it = light_quarks.find(std::abs(pdgId));
        if (it == light_quarks.end()) return false;
        return true;
    };
    p1 = std::find_if(pdgIds_begin, pdgIds_end, IsLightQuark);
    p2 = std::find_if(p1 + 1, pdgIds_end, IsLightQuark);

    idx1 = p1 - pdgIds_begin;
    idx2 = p2 - pdgIds_begin;

    Int_t m1idx, m2idx;
    m1idx = motherIdxs[idx1];
    m2idx = motherIdxs[idx2];

    if (m1idx == m2idx  && std::abs(pdgIds[m1idx]) == 24 && std::abs(pdgIds[m2idx]) == 24 && statuses[idx1] == 23 && statuses[idx2] == 23)
    {
        q1 = idx1;
        q2 = idx2;
    }

    auto IsLepton = [](Int_t pdgId)
    {
        if (std::abs(pdgId) == 11 || std::abs(pdgId) == 13) return true;
        return false;
    };

    auto IsNu = [](Int_t pdgId)
    {
        if (std::abs(pdgId) == 12 || std::abs(pdgId) == 14) return true;
        return false;
    };
    p1 = std::find_if(pdgIds_begin, pdgIds_end, IsLepton);
    p2 = std::find_if(p1 + 1, pdgIds_end, IsNu);

    idx1 = p1 - pdgIds_begin;
    idx2 = p2 - pdgIds_begin;
    m1idx = motherIdxs[idx1];
    m2idx = motherIdxs[idx2];

    if (m1idx == m2idx  && std::abs(pdgIds[m1idx]) == 24 && std::abs(pdgIds[m2idx]) == 24)
    {
        l = idx1;
        nu = idx2;
    }

    return res;
}

GenJetIndex Match(GenPartIndex const& genPartIdx, PtEtaPhiMArray const& genPart, PtEtaPhiMArray const& genJet, Int_t const* genPartId, Int_t const* jetPartFlav)
{
    GenJetIndex res;
    auto& [bj1_idx, bj2_idx, lj1_idx, lj2_idx] = res;
    auto const& [h1, h2, w1, w2, b1, b2, q1, q2, l, nu] = genPartIdx;

    std::vector<Int_t> matched_jets{};
    std::vector<Int_t> target_quarks{b1, b2, q1, q2};

    Int_t nJets = genJet.n;

    for (auto const& tarIdx: target_quarks)
    {
        Int_t pdgId = genPartId[tarIdx];
        std::vector<Int_t> matchable_jets{}; // store indices of jets that can be matched to current quark
        auto IsMatchable = [&pdgId](Int_t part_flav)
        {
            if (std::abs(pdgId) < 5) return (pdgId == part_flav || part_flav == 0 || part_flav == 21);
            return (pdgId == part_flav || part_flav == 0); 
        }; // allow quarks match to jets with undefined flavor (0); can be removed later
        
        auto Insert = [&IsMatchable, &matchable_jets, count = 0](Int_t part_flav) mutable
        {
            if (IsMatchable(part_flav)) matchable_jets.push_back(count);
            count++;
        }; 

        std::for_each(jetPartFlav, jetPartFlav + nJets, Insert); // makes list of indices (pointing to position in GenJet_* arrays) of jets that can be matched

        std::cout << "\tpdgId = " << pdgId << ": ";
        std::copy(matchable_jets.begin(), matchable_jets.end(), std::ostream_iterator<Int_t>(std::cout, " "));
        std::cout << "\n";

        // set target quark momentum
        auto const& [part_pt, part_eta, part_phi, part_m, n_parts] = genPart;
        TLorentzVector part;
        part.SetPtEtaPhiM(part_pt[tarIdx], part_eta[tarIdx], part_phi[tarIdx], part_m[tarIdx]);

        // iterate over mathcable jets to find best match and save its index
        size_t potential_matches = matchable_jets.size();
        Int_t match_idx = -1; // index for the match; if nothing found should stay -1
        Float_t match_dR = 0.4;  
        for (Int_t jetIdx = 0; jetIdx < static_cast<Int_t>(potential_matches); ++jetIdx)
        {
            // set current jet momentum
            auto const& [jet_pt, jet_eta, jet_phi, jet_m, n_jets] = genJet;
            TLorentzVector jet;
            jet.SetPtEtaPhiM(jet_pt[jetIdx], jet_eta[jetIdx], jet_phi[jetIdx], jet_m[jetIdx]);

            Float_t dR = jet.DeltaR(part);
            if (dR > 0.4) 
            {
                std::cout << "\t\tdiscarded:\n";
                std::cout << "\t\t\t";
                Print(part);
                std::cout << "\t\t\t";
                Print(jet);
                std::cout << "\t\t\tdR = " << dR << "\n";
                continue;
            }
            if (dR < match_dR)
            {
                auto it = std::find(matched_jets.begin(), matched_jets.end(), jetIdx);
                if (it != matched_jets.end()) continue;
                match_dR = dR;
                match_idx = jetIdx;
                std::cout << "\t\tattempting to match:\n";
                std::cout << "\t\t\t";
                Print(part);
                std::cout << "\t\t\t";
                Print(jet);
                std::cout << "\t\t\tdR = " << dR << "\n";
            }
            else
            {
                std::cout << "\t\tunable to match:\n";
                std::cout << "\t\t\t";
                Print(part);
                std::cout << "\t\t\t";
                Print(jet);
                std::cout << "\t\t\tdR = " << dR << "\n";
            }   
        }
        matched_jets.push_back(match_idx);
        std::cout << "\tmatches: ";
        std::copy(matched_jets.begin(), matched_jets.end(), std::ostream_iterator<Int_t>(std::cout, " "));
        std::cout << "\n";
    }

    bj1_idx = matched_jets[0];
    bj2_idx = matched_jets[1];
    lj1_idx = matched_jets[2];
    lj2_idx = matched_jets[3];

    return res;
}