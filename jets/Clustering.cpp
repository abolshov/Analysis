#include "Clustering.hpp"

std::vector<Int_t> GetFinalParticles(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart)
{
    std::vector<Int_t> indices(nGenPart);
    std::vector<Bool_t> mask(nGenPart, false);
    for (Int_t pos = 0; pos < static_cast<Int_t>(nGenPart); ++pos)
    {
        Int_t motherIdx = motherIdxs[pos];
        if (motherIdx == -1) continue;
        mask[motherIdx] = true;
    }
    std::generate_n(indices.begin(), nGenPart, [count = 0]() mutable { return count++; });

    auto IsMother = [&mask](Int_t pos) { return mask[pos]; };
    indices.erase(std::remove_if(indices.begin(), indices.end(), IsMother), indices.end());

    return indices;
}

std::vector<Bool_t> PossibleJetConstituents(Int_t const* motherIdxs, Int_t const* StatusArr, Int_t nGenPart, GenPartIndex const& idx)
{
    std::vector<Bool_t> mask(nGenPart, true);
    Int_t lep_pos = idx.l;
    Int_t nu_pos = idx.nu;
    for (Int_t pos = 0; pos < static_cast<Int_t>(nGenPart); ++pos)
    {
        Int_t motherIdx = motherIdxs[pos];
        Int_t status = StatusArr[pos];
        if (IsIncoming(status) || pos == lep_pos || pos == nu_pos) mask[pos] = false;
        if (motherIdx == -1) continue;
        mask[motherIdx] = false;
    }
    return mask;
}