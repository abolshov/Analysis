#ifndef GEN_LVL_CUTS_HPP
#define GEN_LVL_CUTS_HPP

#include "Specification.hpp"
#include "Event.hpp"

class BQuarkAcceptCut final : public Specification<Event>
{
    public:
    BQuarkAcceptCut(Float_t pt, Float_t eta);
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_pt{};
    Float_t m_eta{};
};

class Resolved2bCut final : public Specification<Event>
{
    public:
    explicit Resolved2bCut(Float_t dr_thresh);
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_threshold{};
};


#endif