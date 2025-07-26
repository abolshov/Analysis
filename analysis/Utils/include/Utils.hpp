#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

// implementation taken from https://cpp-optimizations.netlify.app/strings_concatenation/
// with minor correction

//--- functions to calculate the total size ---
size_t StrSize(const char* str) 
{
    return strlen(str);
}

size_t StrSize(const std::string& str) 
{
    return str.size();
}

template <class Head, class... Tail>
size_t StrSize(const Head& head, Tail const&... tail) 
{
    return StrSize(head) + StrSize(tail...);
}

//--- functions to append strings together ---
template <class Head>
void StrAppend(std::string& out, const Head& head) 
{
    out += head;
}

template <class Head, class... Args>
void StrAppend(std::string& out, const Head& head, Args const&... args) 
{
    out += head;
    StrAppend(out, args...);
}

//--- Finally, the function to concatenate strings ---
template <class... Args> 
std::string StrCat(Args const&... args) 
{
    size_t tot_size = StrSize(args...);
    std::string out;
    out.reserve(tot_size);

    StrAppend(out, args...);
    return out;
}

#endif