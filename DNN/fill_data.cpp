#include <iostream>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include <math.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TSpline.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TStyle.h"
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <functional>
#include "TLegend.h"
#include <numeric>
#include "TRandom3.h"
#include <cmath>
#include <regex>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>

float min_dR(TLorentzVector const& part, std::vector<TLorentzVector> const& others)
{
	std::vector<float> tmp_v;
	for (auto const& other: others)
	{
		tmp_v.push_back(part.DeltaR(other));
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_dR(std::vector<TLorentzVector> const& others)
{
	std::vector<float> tmp_v;
	for (size_t i = 0; i < others.size(); ++i)
	{
		for (size_t j = i + 1; j < others.size(); ++j)
		{
			float tmp = others[i].DeltaR(others[j]);
			tmp_v.push_back(tmp);
		}
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_dPhi(TLorentzVector const& part, std::vector<TLorentzVector> const& others)
{
	std::vector<float> tmp_v;
	for (auto const& other: others)
	{
		tmp_v.push_back(part.DeltaPhi(other));
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_dPhi(std::vector<TLorentzVector> const& others)
{
	std::vector<float> tmp_v;
	for (size_t i = 0; i < others.size(); ++i)
	{
		for (size_t j = i + 1; j < others.size(); ++j)
		{
			tmp_v.push_back(others[i].DeltaPhi(others[j]));
		}
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_func(std::vector<TLorentzVector> const& parts, std::function<float(TLorentzVector const& p1, TLorentzVector const& p2)> func)
{
	std::vector<float> tmp_v;
	for (size_t i = 0; i < parts.size(); ++i)
	{
		for (size_t j = 0; j < parts.size(); ++j)
		{
			tmp_v.push_back(func(parts[i], parts[j]));
		}
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_func(TLorentzVector const& part, std::vector<TLorentzVector> const& others, std::function<float(TLorentzVector const& p1, TLorentzVector const& p2)> func)
{
	std::vector<float> tmp_v;
	for (auto const& other: others)
	{
		tmp_v.push_back(func(part, other));
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

std::pair<TLorentzVector, TLorentzVector> select_bjets(std::multimap<float, TLorentzVector> const& jet_data)
{
	// assert (!jet_data.empty() && jet_data.size() != 1);
	if (jet_data.empty() || jet_data.size() == 1)
	{
		throw std::invalid_argument("less than two ak4 jets in the event");
	}

	float const lower_bound_score = 0.3093;
	float const upper_bound_score = 0.7221;

	TLorentzVector jet1, jet2;
	if (jet_data.size() == 2)
	{
		auto jet_it = jet_data.begin();
		jet1 = jet_it->second;
		jet2 = std::next(jet_it)->second;
	}
	else
	{
		std::multimap<float, TLorentzVector>::const_iterator itlow = jet_data.lower_bound(lower_bound_score); 
  		std::multimap<float, TLorentzVector>::const_iterator itup = jet_data.upper_bound(upper_bound_score);

		float min_mass_diff = std::numeric_limits<float>::max();
		if (itlow != jet_data.cend())
		{
			for (auto med_btag_it = itlow; med_btag_it != jet_data.cend(); ++med_btag_it)
			{
				for (auto second_jet_it = jet_data.cbegin(); second_jet_it != jet_data.cend(); ++second_jet_it)
				{
					float cur_mass_diff = abs((med_btag_it->second + second_jet_it->second).M() - 125.0);
					if (cur_mass_diff < min_mass_diff)
					{
						jet1 = med_btag_it->second;
						jet2 = second_jet_it->second;
						min_mass_diff = cur_mass_diff;
					}
				}
			}
		}
		else
		{
			for (auto jet1_it = jet_data.cbegin(); jet1_it != jet_data.cend(); ++jet1_it)
			{
				for (auto jet2_it = std::next(jet1_it); jet2_it != jet_data.cend(); ++jet2_it)
				{
					float cur_mass_diff = abs((jet1_it->second + jet2_it->second).M() - 125.0);
					if (cur_mass_diff < min_mass_diff)
					{
						jet1 = jet1_it->second;
						jet2 = jet2_it->second;
						min_mass_diff = cur_mass_diff;
					}
				}
			}
		}
	}
	return std::make_pair(jet1, jet2);
}

int main() {
	TFile *myFile = TFile::Open("out_radion_2L_2016_m400.root");
	// TFile *myFile = TFile::Open("out_ttbar_2L_2016.root");
	TDirectory *dir = (TDirectory*)myFile;
	TTree *myTree = (TTree*)dir->Get("Double_Tree");

	float ak4_jet0_E;
	myTree->SetBranchAddress("ak4_jet0_E", &ak4_jet0_E);

	float ak4_jet0_phi;
	myTree->SetBranchAddress("ak4_jet0_phi", &ak4_jet0_phi);

	float ak4_jet0_pt;
	myTree->SetBranchAddress("ak4_jet0_pt", &ak4_jet0_pt);

	float ak4_jet0_px;
	myTree->SetBranchAddress("ak4_jet0_px", &ak4_jet0_px);

	float ak4_jet0_py;
	myTree->SetBranchAddress("ak4_jet0_py", &ak4_jet0_py);

	float ak4_jet0_pz;
	myTree->SetBranchAddress("ak4_jet0_pz", &ak4_jet0_pz);

	float ak4_jet1_E;
	myTree->SetBranchAddress("ak4_jet1_E", &ak4_jet1_E);

	float ak4_jet1_phi;
	myTree->SetBranchAddress("ak4_jet1_phi", &ak4_jet1_phi);

	float ak4_jet1_pt;
	myTree->SetBranchAddress("ak4_jet1_pt", &ak4_jet1_pt);

	float ak4_jet1_px;
	myTree->SetBranchAddress("ak4_jet1_px", &ak4_jet1_px);

	float ak4_jet1_py;
	myTree->SetBranchAddress("ak4_jet1_py", &ak4_jet1_py);

	float ak4_jet1_pz;
	myTree->SetBranchAddress("ak4_jet1_pz", &ak4_jet1_pz);

	float ak4_jet2_E;
	myTree->SetBranchAddress("ak4_jet2_E", &ak4_jet2_E);

	float ak4_jet2_phi;
	myTree->SetBranchAddress("ak4_jet2_phi", &ak4_jet2_phi);

	float ak4_jet2_pt;
	myTree->SetBranchAddress("ak4_jet2_pt", &ak4_jet2_pt);

	float ak4_jet2_px;
	myTree->SetBranchAddress("ak4_jet2_px", &ak4_jet2_px);

	float ak4_jet2_py;
	myTree->SetBranchAddress("ak4_jet2_py", &ak4_jet2_py);

	float ak4_jet2_pz;
	myTree->SetBranchAddress("ak4_jet2_pz", &ak4_jet2_pz);

	float ak4_jet3_E;
	myTree->SetBranchAddress("ak4_jet3_E", &ak4_jet3_E);

	float ak4_jet3_phi;
	myTree->SetBranchAddress("ak4_jet3_phi", &ak4_jet3_phi);

	float ak4_jet3_pt;
	myTree->SetBranchAddress("ak4_jet3_pt", &ak4_jet3_pt);

	float ak4_jet3_px;
	myTree->SetBranchAddress("ak4_jet3_px", &ak4_jet3_px);

	float ak4_jet3_py;
	myTree->SetBranchAddress("ak4_jet3_py", &ak4_jet3_py);

	float ak4_jet3_pz;
	myTree->SetBranchAddress("ak4_jet3_pz", &ak4_jet3_pz);

	float ak8_jet0_E;
	myTree->SetBranchAddress("ak8_jet0_E", &ak8_jet0_E);

	float ak8_jet0_phi;
	myTree->SetBranchAddress("ak8_jet0_phi", &ak8_jet0_phi);

	float ak8_jet0_pt;
	myTree->SetBranchAddress("ak8_jet0_pt", &ak8_jet0_pt);

	float ak8_jet0_px;
	myTree->SetBranchAddress("ak8_jet0_px", &ak8_jet0_px);

	float ak8_jet0_py;
	myTree->SetBranchAddress("ak8_jet0_py", &ak8_jet0_py);

	float ak8_jet0_pz;
	myTree->SetBranchAddress("ak8_jet0_pz", &ak8_jet0_pz);

	float lep0_E;
	myTree->SetBranchAddress("lep0_E", &lep0_E);

	int lep0_pdgId;
	myTree->SetBranchAddress("lep0_pdgId", &lep0_pdgId);

	float lep0_phi;
	myTree->SetBranchAddress("lep0_phi", &lep0_phi);

	float lep0_pt;
	myTree->SetBranchAddress("lep0_pt", &lep0_pt);

	float lep0_px;
	myTree->SetBranchAddress("lep0_px", &lep0_px);

	float lep0_py;
	myTree->SetBranchAddress("lep0_py", &lep0_py);

	float lep0_pz;
	myTree->SetBranchAddress("lep0_pz", &lep0_pz);

	float lep1_E;
	myTree->SetBranchAddress("lep1_E", &lep1_E);

	int lep1_pdgId;
	myTree->SetBranchAddress("lep1_pdgId", &lep1_pdgId);

	float lep1_phi;
	myTree->SetBranchAddress("lep1_phi", &lep1_phi);

	float lep1_pt;
	myTree->SetBranchAddress("lep1_pt", &lep1_pt);

	float lep1_px;
	myTree->SetBranchAddress("lep1_px", &lep1_px);

	float lep1_py;
	myTree->SetBranchAddress("lep1_py", &lep1_py);

	float lep1_pz;
	myTree->SetBranchAddress("lep1_pz", &lep1_pz);

	float met_E;
	myTree->SetBranchAddress("met_E", &met_E);

	float met_px;
	myTree->SetBranchAddress("met_px", &met_px);

	float met_py;
	myTree->SetBranchAddress("met_py", &met_py);

	float met_pz;
	myTree->SetBranchAddress("met_pz", &met_pz);

	int dnn_truth_value;
	myTree->SetBranchAddress("dnn_truth_value", &dnn_truth_value);

	float ak4_jet0_btagDeepFlavB;
	myTree->SetBranchAddress("ak4_jet0_btagDeepFlavB", &ak4_jet0_btagDeepFlavB);

	float ak4_jet1_btagDeepFlavB;
	myTree->SetBranchAddress("ak4_jet1_btagDeepFlavB", &ak4_jet1_btagDeepFlavB);

	float ak4_jet2_btagDeepFlavB;
	myTree->SetBranchAddress("ak4_jet2_btagDeepFlavB", &ak4_jet2_btagDeepFlavB);

	float ak4_jet3_btagDeepFlavB;
	myTree->SetBranchAddress("ak4_jet3_btagDeepFlavB", &ak4_jet3_btagDeepFlavB);


	int nEvents = myTree->GetEntries();
    auto branches = myTree->GetListOfBranches();
    int nBranches = myTree->GetNbranches();
	std::vector<std::string> branch_names;
	
	std::regex pat_jet("_jet[0-9]_[pE]");
    std::regex pat_met("met_[pE]");
    std::regex pat_lep("lep[0-9]_[pE]");

	for (size_t i = 0; i < nBranches; ++i)
	{
		std::string bName = branches->At(i)->GetName();
		if (std::regex_search(bName, pat_jet) || 
			std::regex_search(bName, pat_lep) || 
			std::regex_search(bName, pat_met))
		{
			branch_names.push_back(bName);
		}
	}

	std::cout << "# branches = " << branch_names.size() << std::endl;

	// std::ofstream data;
	// data.open("data_signal.csv");
	// data.open("data_backgr.csv");
	// for (auto const& bName: branch_names)
	// {
	// 	data << bName << ",";
	// }
	// data << "ll_dR,ll_dPhi,ll_mass,llmet_dPhi,leadl_ak4jet_mindR,subl_ak4jet_mindR,min_jet_dR,min_jet_dPhi,lmet_mass,";
	// data << "dnn_truth_value" << "\n";

	std::map<int, int> btag_counts;
	btag_counts[0] = 0;
	btag_counts[1] = 0;
	btag_counts[2] = 0;
	btag_counts[3] = 0;
	btag_counts[4] = 0;

	std::map<int, int> zero_momentum_count;
	zero_momentum_count[0] = 0;
	zero_momentum_count[1] = 0;
	zero_momentum_count[2] = 0;
	zero_momentum_count[3] = 0;
	zero_momentum_count[4] = 0;

	int bad_jet_event_count = 0;
	int invalid_jet_event_count = 0;

	TH1F* h_mass = new TH1F("h_mass", "higgs from ak4 jets", 101, 0.0, 250.0);

	for (size_t i = 0; i < nEvents; ++i)
	{
		// i = 298300;
		// myTree->GetEntry(246019);
		myTree->GetEntry(i);

		TLorentzVector l0, l1, met, ak4jet0, ak4jet1, ak4jet2, ak4jet3;
		l0.SetPxPyPzE(lep0_px, lep0_py, lep0_pz, lep0_E);
		l1.SetPxPyPzE(lep1_px, lep1_py, lep1_pz, lep1_E);
		met.SetPxPyPzE(met_px, met_py, met_pz, met_E);
		ak4jet0.SetPxPyPzE(ak4_jet0_px, ak4_jet0_py, ak4_jet0_pz, ak4_jet0_E);
		ak4jet1.SetPxPyPzE(ak4_jet1_px, ak4_jet1_py, ak4_jet1_pz, ak4_jet1_E);
		ak4jet2.SetPxPyPzE(ak4_jet2_px, ak4_jet2_py, ak4_jet2_pz, ak4_jet2_E);
		ak4jet3.SetPxPyPzE(ak4_jet3_px, ak4_jet3_py, ak4_jet3_pz, ak4_jet3_E);

		std::vector<TLorentzVector> jets = {ak4jet0, ak4jet1, ak4jet2, ak4jet3};
		std::vector<TLorentzVector> leptons = {l0, l1};
		std::vector<float> btagScores = {ak4_jet0_btagDeepFlavB, ak4_jet1_btagDeepFlavB, ak4_jet2_btagDeepFlavB, ak4_jet3_btagDeepFlavB};
		std::multimap<float, TLorentzVector> jet_data;

		int n_zero_btags = std::count(btagScores.begin(), btagScores.end(), 0.0);
		btag_counts[4 - n_zero_btags] += 1;

		int n_zero_jets = std::count_if(jets.begin(), jets.end(), [](TLorentzVector const& p) { return p == TLorentzVector(0,0,0,0); });
		zero_momentum_count[n_zero_jets] += 1;

		// int n_zero_lep = std::count_if(leptons.begin(), leptons.end(), [](TLorentzVector const& p) { return p == TLorentzVector(0,0,0,0); });
		// zero_momentum_count[n_zero_lep] += 1;

		bool valid_event = (n_zero_jets < 3 && n_zero_btags < 3);

		if (!valid_event)
		{
			continue;
		}

		// auto print_4vector = [](TLorentzVector const& p) 
		// { 
		// 	std::cout << "(" << p.E() << ", " 
		// 			  << p.Px() << ", " 
		// 			  << p.Py() << ", "
		// 			  << p.Pz() << ")\n"; 
		// };

		// for (size_t j = 0; j < jets.size(); ++j)
		// {
		// 	std::cout << btagScores[j] << ": ";
		// 	print_4vector(jets[j]);
		// }
		// std::cout << "inv_mass = " << (jets[0] + jets[1]).M() << std::endl;
		// std::cout << "-------------------------------------------------\n";

		std::transform(btagScores.begin(), 
					   btagScores.end(), 
					   jets.begin(), 
					   std::inserter(jet_data, jet_data.end()), 
					   [](float a, TLorentzVector b) { return std::make_pair(a, b); });

		jet_data.erase(0); // now it contains only jets with nonzero momentum and btagging score
		
		auto print_jet_data = [](std::pair<float, TLorentzVector> const& p) 
		{ 
			std::cout << p.first << ": (" << p.second.E() << ", " 
					  					  << p.second.Px() << ", " 
										  << p.second.Py() << ", "
										  << p.second.Pz() << ")\n"; 
		};

		// std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
		// std::cout << "-------------------------------------------------\n";
		// break;

		std::pair<TLorentzVector, TLorentzVector> jets_from_h;
		try 
		{
			// std::cout << "Event #" << i << ": jet_data.size() = " << jet_data.size() << std::endl; 
			// std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
			jets_from_h = select_bjets(jet_data);
			float dijet_mass = (jets_from_h.first + jets_from_h.second).M();
			if (dijet_mass < 90.0 || dijet_mass > 160.0)
			{
				++bad_jet_event_count;
				// std::cout << "Event #" << i << ": jet_data.size() = " << jet_data.size() << std::endl; 
				// std::cout << "selected dijet_mass = " << dijet_mass << std::endl;
				// std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
				// std::cout << "-------------------------------------------------\n";
			}
			h_mass->Fill(dijet_mass);
			// std::cout << "-------------------------------------------------\n";
		}
		catch (std::invalid_argument& e)
		{
			valid_event = false;
			std::cout << "Event #" << i << ": " << e.what() << std::endl;
			std::cout << "jet_data.size() = " << jet_data.size() << std::endl;
			std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
			std::cout << "-------------------------------------------------\n";
			++invalid_jet_event_count;
		}

		TLorentzVector jet1, jet2;
		if (valid_event)
		{
			jet1 = jets_from_h.first;
			jet2 = jets_from_h.second;
		}

		float ll_dR = l0.DeltaR(l1);
		float ll_dPhi = l0.DeltaPhi(l1);
		float ll_mass = (l0 + l1).M();
		float llmet_dPhi = met.DeltaPhi(l0 + l1);
		
		TLorentzVector lead_l, sub_l;
		if (l0.Pt() > l1.Pt())
		{
			lead_l = l0;
			sub_l = l1;
		}
		else
		{
			lead_l = l1;
			sub_l = l0;
		}

		auto dR = [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaR(p2); };
		auto dPhi = [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaPhi(p2); };

		// float leadl_ak4jet_mindR = min_dR(lead_l, jets);
		// float subl_ak4jet_mindR = min_dR(sub_l, jets);
		float leadl_ak4jet_mindR = min_func(lead_l, jets, dR);
		float subl_ak4jet_mindR = min_func(sub_l, jets, dR);

		float min_jet_dR = min_dR(jets);
		float min_jet_dPhi = min_dPhi(jets);

		float lmet_mass = (lead_l + met).M();

		// std::cout << "ll_dR = " << ll_dR << "\n" 
		// 		  << "ll_dPhi = " << ll_dPhi << "\n"
		// 		  << "ll_mass = " << ll_mass << "\n"
		// 		  << "llmet_dPhi = " << llmet_dPhi << "\n"
		// 		  << "leadl_ak4jet_mindR = " << leadl_ak4jet_mindR << "\n"
		// 		  << "subl_ak4jet_mindR = " << subl_ak4jet_mindR << "\n"
		// 		  << "min_jet_dR = " << min_jet_dR << "\n"
		// 		  << "min_jet_dPhi = " << min_jet_dPhi << "\n"
		// 		  << "lmet_mass = " << lmet_mass << "\n"
		// 		  << "--------------------------------\n";
	
		// break;

		// data << ak4_jet0_E << ","
		// 	 << ak4_jet0_phi << ","
		// 	 << ak4_jet0_pt << ","
		// 	 << ak4_jet0_px << ","
		// 	 << ak4_jet0_py << ","
		// 	 << ak4_jet0_pz << ","
		// 	 << ak4_jet1_E << ","
		// 	 << ak4_jet1_phi << ","
		// 	 << ak4_jet1_pt << ","
		// 	 << ak4_jet1_px << ","
		// 	 << ak4_jet1_py << ","
		// 	 << ak4_jet1_pz << ","
		// 	 << ak4_jet2_E << ","
		// 	 << ak4_jet2_phi << ","
		// 	 << ak4_jet2_pt << ","
		// 	 << ak4_jet2_px << ","
		// 	 << ak4_jet2_py << ","
		// 	 << ak4_jet2_pz << ","
		// 	 << ak4_jet3_E << ","
		// 	 << ak4_jet3_phi << ","
		// 	 << ak4_jet3_pt << ","
		// 	 << ak4_jet3_px << ","
		// 	 << ak4_jet3_py << ","
		// 	 << ak4_jet3_pz << ","
		// 	 << ak8_jet0_E << ","
		// 	 << ak8_jet0_phi << ","
		// 	 << ak8_jet0_pt << ","
		// 	 << ak8_jet0_px << ","
		// 	 << ak8_jet0_py << ","
		// 	 << ak8_jet0_pz << ","
		// 	 << lep0_E << ","
		// 	 << lep0_pdgId << ","
		// 	 << lep0_phi << ","
		// 	 << lep0_pt << ","
		// 	 << lep0_px << ","
		// 	 << lep0_py << ","
		// 	 << lep0_pz << ","
		// 	 << lep1_E << ","
		// 	 << lep1_pdgId << ","
		// 	 << lep1_phi << ","
		// 	 << lep1_pt << ","
		// 	 << lep1_px << ","
		// 	 << lep1_py << ","
		// 	 << lep1_pz << ","
		// 	 << met_E << ","
		// 	 << met_px << ","
		// 	 << met_py << ","
		// 	 << met_pz << ","
		//   << ll_dR << ","
		// 	 << ll_dPhi << ","
		// 	 << ll_mass << ","
		// 	 << llmet_dPhi << ","
		// 	 << leadl_ak4jet_mindR << ","
		// 	 << subl_ak4jet_mindR << ","
		// 	 << min_jet_dR << ","
		// 	 << min_jet_dPhi << ","
		// 	 << lmet_mass << ","
		// 	 << dnn_truth_value << "\n";
		// break;
	}

	TCanvas* c1 = new TCanvas("c1", "c1");
    c1->SetCanvasSize(1500, 1500);
    c1->SetWindowSize(1000, 1000);

	h_mass->SetLineWidth(3);
	h_mass->DrawNormalized();
	c1->SaveAs("h_mass.png");

	std::for_each(btag_counts.begin(), btag_counts.end(), [](std::pair<int, int> const& p) { std::cout << p.first << " b-tagged jet(s): " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::for_each(zero_momentum_count.begin(), zero_momentum_count.end(), [](std::pair<int, int> const& p) { std::cout << p.first << " zero jet(s): " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::cout << "b jets not selected: invalid_jet_event_count = " << invalid_jet_event_count << std::endl;
	std::cout << "higgs from bb mass wrong: bad_jet_event_count = " << bad_jet_event_count << std::endl;
	std::cout << "---------------------------\n";
	// data.close();
	return 0;
}
