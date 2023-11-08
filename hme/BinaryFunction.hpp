#ifndef BIN_FUNC_HPP
#define BIN_FUNC_HPP

#include <functional>
#include <vector>

#include "TLorentzVector.h"

template <typename Comparator = std::less<>>
struct BinaryFunction
{
    using Function = std::function<float(TLorentzVector const& v1, TLorentzVector const& v2)>;

    BinaryFunction(std::string const& name, Function func, float cut, Comparator cmp = Comparator{}) noexcept
    : m_name(name), m_func(func), m_cut(cut), m_cmp(cmp), m_counter(0) {}

    ~BinaryFunction() = default;

    BinaryFunction(BinaryFunction const& other) = delete;
    BinaryFunction& operator=(BinaryFunction const& other) = delete;
    BinaryFunction(BinaryFunction&& other) = delete;
    BinaryFunction& operator=(BinaryFunction&& other) = delete;

    inline float operator()(TLorentzVector const& v1, TLorentzVector const& v2) const { return m_func(v1, v2); }
    inline std::string const& Name() const { return m_name; }
    inline bool IsValid(TLorentzVector const& v1, TLorentzVector const& v2) const { return m_cmp(m_func(v1, v2), m_cut); }

    void Print(std::vector<TLorentzVector> const& particles)
    {
        std::cout << "Printing " << m_name << ":\n";
        size_t bad_count = 0;
        for (size_t i = 0; i < particles.size(); ++i)
        {
            for (size_t j = i + 1; j < particles.size(); ++j)
            {
                if (IsValid(particles[i], particles[j])) 
                {
                    std::cout << "\t" << m_name << "(" << i << ", " << j << "): " << m_func(particles[i], particles[j]) << "\n";
                    ++bad_count;
                }
            }
        }

        if (bad_count > 4) 
        {
            std::cout << "WARNING: found " << bad_count << " pairs not passing the cut\n";
            ++m_counter;
        }
    }

    inline size_t Count() const { return m_counter; }

    private:
    std::string m_name;
    Function m_func;
    float m_cut;
    Comparator m_cmp;
    size_t m_counter;
};

#endif