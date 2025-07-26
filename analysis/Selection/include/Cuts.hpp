#ifndef CUTS_HPP
#define CUTS_HPP

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

#endif