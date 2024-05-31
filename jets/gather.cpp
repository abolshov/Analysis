#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "MatchingTools.hpp"

static constexpr int MAX_AK4_GENJET = 21;
static constexpr int MAX_GENPART = 270;

int main([[maybe_unused]] int argc, char* argv[])
{
    auto start = std::chrono::system_clock::now();
    std::string input_file_name = argv[1];
    std::cout << "Processing " << input_file_name << "\n";
    std::unique_ptr<TFile> input_file(TFile::Open(input_file_name.c_str()));

    auto it = std::find(input_file_name.begin(), input_file_name.end(), '_');
    std::string output_file_name("JetNetTrain");
    output_file_name.append(it, input_file_name.end());
    std::cout << "Will write to " << output_file_name << "\n";
    std::unique_ptr<TFile> output_file(TFile::Open(output_file_name.c_str(), "RECREATE"));

    auto input_tree = std::unique_ptr<TTree>(static_cast<TTree*>(input_file->Get("Events")));
    auto output_tree = std::make_unique<TTree>("JetNetTree", "Tree containing data for JetNet");

    UInt_t          nGenPart;
    Float_t         GenPart_eta[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_mass[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_phi[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_pt[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_pdgId[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_status[MAX_GENPART];   //[nGenPart]

    // ak4 jets
    UInt_t          nGenJetAK4;
    Float_t         GenJetAK4_eta[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_mass[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_phi[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_pt[MAX_AK4_GENJET];   //[nGenJetAK4]
    Int_t           GenJetAK4_partonFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    UChar_t         GenJetAK4_hadronFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]

    Float_t         GenMET_phi;
    Float_t         GenMET_pt;


    // gen particles 
    input_tree->SetBranchAddress("nGenPart", &nGenPart);
    input_tree->SetBranchAddress("GenPart_eta", &GenPart_eta);
    input_tree->SetBranchAddress("GenPart_mass", &GenPart_mass);
    input_tree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    input_tree->SetBranchAddress("GenPart_pt", &GenPart_pt);
    input_tree->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    input_tree->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    input_tree->SetBranchAddress("GenPart_status", &GenPart_status);
    input_tree->SetBranchAddress("GenPart_phi", &GenPart_phi);

    // ak4 gen jets
    input_tree->SetBranchAddress("nGenJet", &nGenJetAK4);
    input_tree->SetBranchAddress("GenJet_eta", &GenJetAK4_eta);
    input_tree->SetBranchAddress("GenJet_mass", &GenJetAK4_mass);
    input_tree->SetBranchAddress("GenJet_phi", &GenJetAK4_phi);
    input_tree->SetBranchAddress("GenJet_pt", &GenJetAK4_pt);
    input_tree->SetBranchAddress("GenJet_partonFlavour", &GenJetAK4_partonFlavour);
    input_tree->SetBranchAddress("GenJet_hadronFlavour", &GenJetAK4_hadronFlavour);

    input_tree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    input_tree->SetBranchAddress("GenMET_pt", &GenMET_pt);

    Float_t genbjet1_px;
    Float_t genbjet1_py;
    Float_t genbjet1_pz;
    Float_t genbjet1_E;

    Float_t genbjet2_px;
    Float_t genbjet2_py;
    Float_t genbjet2_pz;
    Float_t genbjet2_E;

    Float_t X_px;
    Float_t X_py;
    Float_t X_pz;
    Float_t X_E;

    Float_t X_mass;

    Float_t H_bb_px;
    Float_t H_bb_py;
    Float_t H_bb_pz;
    Float_t H_bb_E;

    Float_t b_px;
    Float_t b_py;
    Float_t b_pz;
    Float_t b_E;

    Float_t bbar_px;
    Float_t bbar_py;
    Float_t bbar_pz;
    Float_t bbar_E;

    Float_t H_WW_px;
    Float_t H_WW_py;
    Float_t H_WW_pz;
    Float_t H_WW_E;

    Float_t LepWfirst_px;
    Float_t LepWfirst_py;
    Float_t LepWfirst_pz;
    Float_t LepWfirst_E;

    Float_t LepWlast_px;
    Float_t LepWlast_py;
    Float_t LepWlast_pz;
    Float_t LepWlast_E;

    Float_t l_px;
    Float_t l_py;
    Float_t l_pz;
    Float_t l_E;

    Float_t nu_px;
    Float_t nu_py;
    Float_t nu_pz;
    Float_t nu_E;

    Float_t HadWfirst_px;
    Float_t HadWfirst_py;
    Float_t HadWfirst_pz;
    Float_t HadWfirst_E;

    Float_t HadWlast_px;
    Float_t HadWlast_py;
    Float_t HadWlast_pz;
    Float_t HadWlast_E;

    Float_t q1_px;
    Float_t q1_py;
    Float_t q1_pz;
    Float_t q1_E;

    Float_t q2_px;
    Float_t q2_py;
    Float_t q2_pz;
    Float_t q2_E;

    output_tree->Branch("genbjet1_px", &genbjet1_px);
    output_tree->Branch("genbjet1_py", &genbjet1_py);
    output_tree->Branch("genbjet1_pz", &genbjet1_pz);
    output_tree->Branch("genbjet1_E", &genbjet1_E);

    output_tree->Branch("genbjet2_px", &genbjet2_px);
    output_tree->Branch("genbjet2_py", &genbjet2_py);
    output_tree->Branch("genbjet2_pz", &genbjet2_pz);
    output_tree->Branch("genbjet2_E", &genbjet2_E);

    output_tree->Branch("X_px", &X_px);
    output_tree->Branch("X_py", &X_py);
    output_tree->Branch("X_pz", &X_pz);
    output_tree->Branch("X_E", &X_E);

    output_tree->Branch("X_mass", &X_mass);

    output_tree->Branch("H_bb_px", &H_bb_px);
    output_tree->Branch("H_bb_py", &H_bb_py);
    output_tree->Branch("H_bb_pz", &H_bb_pz);
    output_tree->Branch("H_bb_E", &H_bb_E);

    output_tree->Branch("b_px", &b_px);
    output_tree->Branch("b_py", &b_py);
    output_tree->Branch("b_pz", &b_pz);
    output_tree->Branch("b_E", &b_E);

    output_tree->Branch("bbar_px", &bbar_px);
    output_tree->Branch("bbar_py", &bbar_py);
    output_tree->Branch("bbar_pz", &bbar_pz);
    output_tree->Branch("bbar_E", &bbar_E);

    output_tree->Branch("H_WW_px", &H_WW_px);
    output_tree->Branch("H_WW_py", &H_WW_py);
    output_tree->Branch("H_WW_pz", &H_WW_pz);
    output_tree->Branch("H_WW_E", &H_WW_E);

    output_tree->Branch("LepWfirst_px", &LepWfirst_px);
    output_tree->Branch("LepWfirst_py", &LepWfirst_py);
    output_tree->Branch("LepWfirst_pz", &LepWfirst_pz);
    output_tree->Branch("LepWfirst_E", &LepWfirst_E);

    output_tree->Branch("LepWlast_px", &LepWlast_px);
    output_tree->Branch("LepWlast_py", &LepWlast_py);
    output_tree->Branch("LepWlast_pz", &LepWlast_pz);
    output_tree->Branch("LepWlast_E", &LepWlast_E);

    output_tree->Branch("l_px", &l_px);
    output_tree->Branch("l_py", &l_py);
    output_tree->Branch("l_pz", &l_pz);
    output_tree->Branch("l_E", &l_E);

    output_tree->Branch("nu_px", &nu_px);
    output_tree->Branch("nu_py", &nu_py);
    output_tree->Branch("nu_pz", &nu_pz);
    output_tree->Branch("nu_E", &nu_E);

    output_tree->Branch("HadWfirst_px", &HadWfirst_px);
    output_tree->Branch("HadWfirst_py", &HadWfirst_py);
    output_tree->Branch("HadWfirst_pz", &HadWfirst_pz);
    output_tree->Branch("HadWfirst_E", &HadWfirst_E);

    output_tree->Branch("HadWlast_px", &HadWlast_px);
    output_tree->Branch("HadWlast_py", &HadWlast_py);
    output_tree->Branch("HadWlast_pz", &HadWlast_pz);
    output_tree->Branch("HadWlast_E", &HadWlast_E);

    output_tree->Branch("q1_px", &q1_px);
    output_tree->Branch("q1_py", &q1_py);
    output_tree->Branch("q1_pz", &q1_pz);
    output_tree->Branch("q1_E", &q1_E);

    output_tree->Branch("q2_px", &q2_px);
    output_tree->Branch("q2_py", &q2_py);
    output_tree->Branch("q2_pz", &q2_pz);
    output_tree->Branch("q2_E", &q2_E);

    int nEvents = input_tree->GetEntries();

    std::cout << std::boolalpha;
    int true_X_mass = 0; 

    for (int i = 0; i < nEvents; ++i)
    {
        input_tree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            true_X_mass = static_cast<int>(GenPart_mass[sig[SIG::X]]);

            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                // CUTS START
                if (!PassLeptonCut(l_p4))
                {
                    continue;
                }

                std::vector<int> selected_jets = PrimaryJetSelection(genjet_ak4);
                if (selected_jets.size() < 4)
                {
                    continue;
                }

                // perform jet cleaning: reject jets overlapping with the lepton
                auto JetLepOverlap = [&l_p4, &genjet_ak4](int i)
                { 
                    TLorentzVector p = GetP4(genjet_ak4, i);
                    return p.DeltaR(l_p4) < DR_THRESH;
                };
                selected_jets.erase(std::remove_if(selected_jets.begin(), selected_jets.end(), JetLepOverlap), selected_jets.end());
                if (selected_jets.size() < 4)
                {
                    continue;
                }

                std::vector<TLorentzVector> b_jets;
                for (auto j: selected_jets)
                {
                    unsigned flav = static_cast<unsigned>(GenJetAK4_hadronFlavour[j]);
                    if (flav == 5)
                    {
                        b_jets.push_back(GetP4(genjet_ak4, j));
                    }
                }

                int num_bflav_jets = b_jets.size();
                if (num_bflav_jets < 1)
                {
                    continue;
                }

                if (num_bflav_jets < 2)
                {
                    continue;
                }

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_B_GENJET_PT; };
                int num_bjets_passing_pt = std::count_if(b_jets.begin(), b_jets.end(), PtCut);
                if (num_bjets_passing_pt < 2)
                {
                    continue;
                }

                auto EtaCut = [](TLorentzVector const& p) { return std::abs(p.Eta()) < MAX_GENJET_ETA; };
                int num_bjets_passing_eta = std::count_if(b_jets.begin(), b_jets.end(), EtaCut);
                if (num_bjets_passing_eta < 2)
                {
                    continue;
                }

                TLorentzVector lq1_p4 = GetP4(genpart, sig[SIG::q1]);
                TLorentzVector lq2_p4 = GetP4(genpart, sig[SIG::q2]);
                TLorentzVector bq_p4 = GetP4(genpart, sig[SIG::b]);
                TLorentzVector bbarq_p4 = GetP4(genpart, sig[SIG::bbar]);

                TLorentzVector X_p4 = GetP4(genpart, sig[SIG::X]);
                TLorentzVector H_bb_p4 = GetP4(genpart, sig[SIG::H_bb]);
                TLorentzVector H_WW_p4 = GetP4(genpart, sig[SIG::H_WW]);

                TLorentzVector LepWfirst_p4 = GetP4(genpart, sig[SIG::LepWfirst]);
                TLorentzVector LepWlast_p4 = GetP4(genpart, sig[SIG::LepWlast]);
                TLorentzVector HadWfirst_p4 = GetP4(genpart, sig[SIG::HadWfirst]);
                TLorentzVector HadWlast_p4 = GetP4(genpart, sig[SIG::HadWlast]);

                TLorentzVector bj1_p4, bj2_p4;
                if (b_jets.size() == 2)
                {
                    bj1_p4 = b_jets[0];
                    bj2_p4 = b_jets[1];
                }
                else 
                {
                    continue;
                }

                genbjet1_px = bj1_p4.Px();
                genbjet1_py = bj1_p4.Py();
                genbjet1_pz = bj1_p4.Pz();
                genbjet1_E = bj1_p4.E();

                genbjet2_px = bj2_p4.Px();
                genbjet2_py = bj2_p4.Py();
                genbjet2_pz = bj2_p4.Pz();
                genbjet2_E = bj2_p4.E();

                X_px = X_p4.Px();
                X_py = X_p4.Py();
                X_pz = X_p4.Pz();
                X_E = X_p4.E();

                X_mass = X_p4.M();

                H_bb_px = H_bb_p4.Px();
                H_bb_py = H_bb_p4.Py();
                H_bb_pz = H_bb_p4.Pz();
                H_bb_E = H_bb_p4.E();

                b_px = bq_p4.Px();
                b_py = bq_p4.Py();
                b_pz = bq_p4.Pz();
                b_E = bq_p4.E();

                bbar_px = bbarq_p4.Px();
                bbar_py = bbarq_p4.Py();
                bbar_pz = bbarq_p4.Pz();
                bbar_E = bbarq_p4.E();

                H_WW_px = H_WW_p4.Px();
                H_WW_py = H_WW_p4.Py();
                H_WW_pz = H_WW_p4.Pz();
                H_WW_E = H_WW_p4.E();

                HadWfirst_px = HadWfirst_p4.Px();
                HadWfirst_py = HadWfirst_p4.Py();
                HadWfirst_pz = HadWfirst_p4.Pz();
                HadWfirst_E = HadWfirst_p4.E();

                HadWlast_px = HadWlast_p4.Px();
                HadWlast_py = HadWlast_p4.Py();
                HadWlast_pz = HadWlast_p4.Pz();
                HadWlast_E = HadWlast_p4.E();

                q1_px = lq1_p4.Px();
                q1_py = lq1_p4.Py();
                q1_pz = lq1_p4.Pz();
                q1_E = lq1_p4.E();

                q2_px = lq2_p4.Px();
                q2_py = lq2_p4.Py();
                q2_pz = lq2_p4.Pz();
                q2_E = lq2_p4.E();

                LepWfirst_px = LepWfirst_p4.Px();
                LepWfirst_py = LepWfirst_p4.Py();
                LepWfirst_pz = LepWfirst_p4.Pz();
                LepWfirst_E = LepWfirst_p4.E();

                LepWlast_px = LepWlast_p4.Px();
                LepWlast_py = LepWlast_p4.Py();
                LepWlast_pz = LepWlast_p4.Pz();
                LepWlast_E = LepWlast_p4.E();

                l_px = l_p4.Px();
                l_py = l_p4.Py();
                l_pz = l_p4.Pz();
                l_E = l_p4.E();

                nu_px = nu_p4.Px();
                nu_py = nu_p4.Py();
                nu_pz = nu_p4.Pz();
                nu_E = nu_p4.E();

                output_tree->Fill();
            }
        }
    }

    output_tree->Write("", TObject::kOverwrite);
    std::cout << "Output tree contains " << output_tree->GetEntries() << " entries\n";

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing sample with mass M = " << true_X_mass << "\n";
    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n";
    return 0;
}