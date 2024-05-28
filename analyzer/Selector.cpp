#include "Selector.hpp"

std::optional<SignalData> Selector::operator()(std::unique_ptr<EventData> const& data)
{
    return std::nullopt;
}