#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include <optional>

#include "SignalData.hpp"

struct Selector 
{
    std::optional<SignalData> operator()(std::unique_ptr<EventData> const& data);

    private:
    std::vector<int> GetSignal(int const* pdg_ids, int const* mothers, int n_gen_part) const;
};

#endif