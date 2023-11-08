#ifndef MYLIB_H
#define MYLIB_H

#include <vector>

namespace bits
{
    size_t extract_bit(int number, size_t bit_id) { return ((number & (1 << (bit_id - 1))) >> (bit_id - 1)); }
}

namespace strings
{
    std::vector<std::string> split(std::string const& s, std::string const& delimiter); 
    std::vector<std::string> split(std::string const& s, char delimeter);
}

#endif