#ifndef TRACER_HPP
#define TRACER_HPP

#include <fstream>

class Tracer
{
    public:
    static Tracer& instance();
    ~Tracer();

    Tracer(Tracer const& other) = delete;
    Tracer(Tracer&& other) = delete;
    Tracer& operator=(Tracer const& other) = delete;
    Tracer& operator=(Tracer&& other) = delete;

    void write(std::string const& msg);

    private:
    std::ofstream m_file;

    Tracer();
};

#endif