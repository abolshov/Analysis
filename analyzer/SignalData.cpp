#include "SignalData.hpp"

explicit SignalData::SignalData(std::unique_ptr<EventData> const& data) 
: m_data(data),
  m_genpart_selection(N_SIG_PART),
  m_bjet_selection(N_B_JETS),
  m_light_jet_selection(N_LIGHT_JETS),
  m_genpart_momenta(),
  m_bjet_momenta(),
  m_light_jet_momenta()  
{}

void SignalData::ComputeMomenta()
{
    if (!m_genpart_selection.empty())
    {
        for (auto idx: m_genpart_selection)
        {
            m_genpart_momenta.push_back(m_data->ComputeP4(idx, m_data->GenPart_pt[idx], m_data->GenPart_eta[idx], m_data->GenPart_phi[idx], m_data->GenPart_mass[idx]));
        }
    }

    if (!m_bjet_selection.empty())
    {
        for (auto idx: m_bjet_selection)
        {
            m_bjet_momenta.push_back(m_data->ComputeP4(idx, m_data->GenJetAK4_pt[idx], m_data->GenJetAK4_eta[idx], m_data->GenJetAK4_phi[idx], m_data->GenJetAK4_mass[idx]));
        }
    }

    if (!m_light_jet_selection.empty())
    {
        for (auto idx: m_light_jet_selection)
        {
            m_light_jet_momenta.push_back(m_data->ComputeP4(idx, m_data->GenJetAK4_pt[idx], m_data->GenJetAK4_eta[idx], m_data->GenJetAK4_phi[idx], m_data->GenJetAK4_mass[idx]));
        }
    }
}

