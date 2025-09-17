#include <iostream>

#include "Event.hpp"
#include "GenLevelCuts.hpp"

int main()
{
    Event event{};
    std::cout << event.n_reco_jet << "\n";
    static constexpr size_t EVT_SZ = sizeof(Event);
    std::cout << "Event size: " << EVT_SZ << "\n";

    BQuarkAcceptCut bq_accept{20.0, 2.4};
    std::cout << bq_accept << " satisfied: " << std::boolalpha << bq_accept.IsSatisfied(event) << "\n";

    ResolvedBQuarksCut res2b{0.8};
    std::cout << res2b << " satisfied: " << std::boolalpha << res2b.IsSatisfied(event) << "\n";

    return 0;
}