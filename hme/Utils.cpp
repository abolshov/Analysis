#include "Utils.hpp"

bool ValidEvent(std::vector<TLorentzVector> const& particles)
{
    auto zeroIt = std::find(particles.begin(), particles.end(), zero);
    if (zeroIt != particles.end()) return false;
    if (particles[PART::b1] == particles[PART::b2]) return false;
    if (particles[PART::q1] == particles[PART::q2]) return false;
    return true;
}

bool PartonCut(std::vector<TLorentzVector> const& particles)
{
    bool good_lq_1 = (particles[PART::q1].Pt() > 20.0f && abs(particles[PART::q1].Eta()) < 2.4f);
    bool good_lq_2 = (particles[PART::q2].Pt() > 20.0f && abs(particles[PART::q2].Eta()) < 2.4f); 
    bool good_light_quarks = (good_lq_1 && good_lq_2);

    bool good_b_1 = (particles[PART::b1].Pt() > 20.0f && abs(particles[PART::b1].Eta()) < 2.4f);
    bool good_b_2 = (particles[PART::b2].Pt() > 20.0f && abs(particles[PART::b2].Eta()) < 2.4f);
    bool good_b_quarks = (good_b_1 && good_b_2);

    bool good_l = (particles[PART::l].Pt() > 10.0f && abs(particles[PART::l].Eta()) < 2.5f);
    return (good_light_quarks && good_b_quarks && good_l);
}

bool GenJetCut(std::vector<TLorentzVector> const& particles)
{
    bool good_lj_1 = (particles[PART::q1].Pt() > 25.0f && abs(particles[PART::q1].Eta()) < 2.4f);
    bool good_lj_2 = (particles[PART::q2].Pt() > 25.0f && abs(particles[PART::q2].Eta()) < 2.4f); 
    bool good_light_jets = (good_lj_1 && good_lj_2);

    bool good_b_1 = (particles[PART::b1].Pt() > 25.0f && abs(particles[PART::b1].Eta()) < 2.4f);
    bool good_b_2 = (particles[PART::b2].Pt() > 25.0f && abs(particles[PART::b2].Eta()) < 2.4f);
    bool good_b_jets = (good_b_1 && good_b_2);

    bool good_l = (particles[PART::l].Pt() > 10.0f && abs(particles[PART::l].Eta()) < 2.5f);
    return (good_light_jets && good_b_jets && good_l);
}

void Print(TLorentzVector const& p, bool EXYZ)
{
    if (EXYZ) 
    {
        std::cout << "(" << p.E() << ", " << p.X() << ", " << p.Y() << ", " << p.Z() << ")\n";
    }
    else
    {
        std::cout << "(" << p.Pt() << ", " << p.Eta() << ", " << p.Phi() << ", " << p.M() << ")\n";
    }
}

bool HasZeroParticle(std::vector<TLorentzVector> const& particles)
{
    TLorentzVector const zero(0.0f, 0.0f, 0.0f, 0.0f);
    auto zeroIt = std::find(particles.begin(), particles.end(), zero);
    if (zeroIt != particles.end()) return true;
    return false;
}

bool IsIdenticalPair(TLorentzVector const& p1, TLorentzVector const& p2)
{
    return p1 == p2;
}
