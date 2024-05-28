#ifndef SIGNAL_DATA_HPP
#define SIGNAL_DATA_HPP

#include <memory>
#include <vector>

#include "TLorentzVector.h"

#include "EventData.hpp"
#include "Constants.hpp"

// class holding information selected by the Selector:
// - indices of signal particles in arrays inside of EventData
// - indices of selected jets
// - 4-momenta of selected partons and jets

class SignalData 
{
    public:
    explicit SignalData(std::unique_ptr<EventData> const& data);

    void ComputeMomenta();
    inline std::vector<int> const& GetSignalIndices() const { return m_genpart_selection; }
    inline std::vector<TLorentzVector>& GetGenpartP4() { return m_genpart_momenta; }
    inline std::vector<TLorentzVector>& GetBJetP4() { return m_bjet_momenta; }
    inline std::vector<TLorentzVector>& GetLightJetP4() { return m_light_jet_momenta; }

    private:
    std::unique_ptr<EventData> const& m_data;

    std::vector<int> m_genpart_selection;
    std::vector<int> m_bjet_selection;
    std::vector<int> m_light_jet_selection;

    std::vector<TLorentzVector> m_genpart_momenta;
    std::vector<TLorentzVector> m_bjet_momenta;
    std::vector<TLorentzVector> m_light_jet_momenta;
};

#endif