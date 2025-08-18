#ifndef GEN_RECO_MATCH_HPP
#define GEN_RECO_MATCH_HPP

#include "Specification.hpp"
#include "Event.hpp"

// inline constexpr Float_t AK4_MATCH_THRSH = 0.4f;
// inline constexpr Float_t AK8_MATCH_THRSH = 0.8f;

class AK4BJetMatchingCut final : public Specification<Event>
{
    public:
    explicit AK4BJetMatchingCut(Float_t dr_thresh);
    ~AK4BJetMatchingCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_match_thresh{};
}; 

class AK8BJetMatchingCut final : public Specification<Event>
{
    public:
    explicit AK8BJetMatchingCut(Float_t dr_thresh);
    ~AK8BJetMatchingCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_match_thresh{};
}; 

#endif