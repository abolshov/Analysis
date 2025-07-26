#include <iostream>

#include "Event.hpp"
#include "Cuts.hpp"

int main()
{
    Event event{};
    std::cout << event.n_reco_jet << "\n";
    static constexpr size_t EVT_SZ = sizeof(Event);
    std::cout << "Event size: " << EVT_SZ << "\n";

    BQuarkAcceptCut cut{20.0, 2.4};
    std::cout << cut << " satisfied: " << std::boolalpha << cut.IsSatisfied(event) << "\n";

    return 0;
}