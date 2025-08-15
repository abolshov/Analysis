#ifndef RECO_LVL_CUTS_HPP
#define RECO_LVL_CUTS_HPP

#include "Specification.hpp"
#include "Event.hpp"

class AK4JetAcceptCut final : public Specification<Event>
{
    AK4JetAcceptCut(Float_t pt, Float_t eta, Int_t num_jets);
    ~AK4JetAcceptCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_pt{};
    Float_t m_eta{};
    Int_t m_num_jets{};
}; 

// right now looks like code duplication, but can be easily extended if btag requirement needs to be included
// ak4 and ak8 jets use different taggers from different branches
class AK8JetAcceptCut final : public Specification<Event>
{
    AK8JetAcceptCut(Float_t pt, Float_t eta, Int_t num_jets);
    ~AK8JetAcceptCut() = default;
    bool IsSatisfied(Event const& event) override;

    private:
    Float_t m_pt{};
    Float_t m_eta{};
    Int_t m_num_jets{};
}; 

#endif