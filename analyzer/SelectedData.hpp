#ifndef SELECTED_DATA_HPP
#define SELECTED_DATA_HPP

#include <memory>
#include <vector>

#include "TLorentzVector.h"

#include "EventData.hpp"
#include "Constants.hpp"

class SelectedData 
{
    public:
    explicit SelectedData(std::unique_ptr<EventData> const& data);

    void ComputeMomenta();
    inline std::vector<int> const& GetSignalIndices() const { return m_genpart_selection; }
    inline std::vector<TLorentzVector>& GetGenpartP4() { return m_genpart_momenta; }
    inline std::vector<TLorentzVector>& GetGenjetP4() { return m_genjet_momenta; }

    private:
    std::unique_ptr<EventData> const& m_data;
    std::vector<int> m_genpart_selection;
    std::vector<int> m_bjet_selection;
    std::vector<int> m_light_jet_selection;
    std::vector<TLorentzVector> m_genpart_momenta;
    std::vector<TLorentzVector> m_bjet_momenta;
};

#endif