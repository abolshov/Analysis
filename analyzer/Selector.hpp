#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include <optional>

#include "SelectedData.hpp"

struct Selector 
{
    std::optional<SelectedData> operator()(std::unique_ptr<EventData> const& data);
};

#endif