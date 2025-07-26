#ifndef SPEC_HPP
#define SPEC_HPP

#include <iostream>

template <typename T>
class Specification
{
    public:
    Specification() = default;

    Specification(std::string name, std::string description)
    :   m_name(std::move(name))
    ,   m_description(std::move(description))
    {}

    virtual bool IsSatisfied(T const& item) = 0;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, Specification<U> const& spec)
    {
        os << spec.m_name << ": " << spec.m_description;
        return os;
    }

    private:
    std::string m_name;
    std::string m_description;
};

#endif