#ifndef PREFETCHER_HPP
#define PREFETCHER_HPP

#include <vector>

#include "Event.hpp"

class Prefetcher
{
    public:
    

    private:
    Event m_buffer{};
    std::vector<Event> m_events;
};

#endif