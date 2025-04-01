#ifndef ESTIM_DATA_HPP
#define ESTIM_DATA_HPP

struct EstimatorData
{
    ULong64_t event_id{};
    Float_t mass{};
    Float_t integral{};
    Float_t width{};
    Float_t peak_value{};
};

#endif