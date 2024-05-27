#ifndef SELECTOR_HPP
#define SELECTOR_HPP

class Selector 
{
    public:
    Selector();
    bool operator()(EventData const& data);
};

#endif