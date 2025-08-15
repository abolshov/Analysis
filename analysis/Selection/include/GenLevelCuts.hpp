#ifndef GEN_LVL_CUTS_HPP
#define GEN_LVL_CUTS_HPP

#include "Specification.hpp"
#include "Event.hpp"

class BQuarkAcceptCut final : public Specification<Event>
{
    public:
    BQuarkAcceptCut(Float_t pt, Float_t eta);
    ~BQuarkAcceptCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_pt{};
    Float_t m_eta{};
};

class Resolved2bCut final : public Specification<Event>
{
    public:
    explicit Resolved2bCut(Float_t dr_thresh);
    ~Resolved2bCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_threshold{};
};

class LeptonAcceptCutDL final : public Specification<Event>
{
    public:
    explicit LeptonAcceptCutDL(Float_t pt);
    ~LeptonAcceptCutDL() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_pt{};
};

#endif