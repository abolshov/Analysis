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


float min_func(std::vector<TLorentzVector> const& parts, std::function<float(TLorentzVector const& p1, TLorentzVector const& p2)> func)
{
	TLorentzVector const zero(0.0, 0.0, 0.0, 0.0);

	for (auto& particle: parts)
	{
		assert(particle != zero);
	}

	assert(parts.size() > 1);

	std::vector<float> tmp_v;
	for (size_t i = 0; i < parts.size(); ++i)
	{
		for (size_t j = i + 1; j < parts.size(); ++j)
		{
			tmp_v.push_back(func(parts[i], parts[j]));
		}
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

float min_func(TLorentzVector const& part, std::vector<TLorentzVector> const& others, std::function<float(TLorentzVector const& p1, TLorentzVector const& p2)> func)
{
	TLorentzVector const zero(0.0, 0.0, 0.0, 0.0);

	assert(part != zero);

	for (auto& particle: others)
	{
		assert(particle != zero);
	}

	assert(others.size() > 1);

	std::vector<float> tmp_v;
	for (auto const& other: others)
	{
		tmp_v.push_back(func(part, other));
	}
	return *std::min_element(tmp_v.begin(),tmp_v.end());
}

std::pair<TLorentzVector, TLorentzVector> select_bjets(std::multimap<float, TLorentzVector> const& jet_data)
{
	
	if (std::count_if(jet_data.begin(), jet_data.end(), 
	    [](std::pair<float, TLorentzVector> const& pair) 
		{ return (pair.second == TLorentzVector(0.0, 0.0, 0.0, 0.0) || pair.first == 0.0); }) != 0)
	{
		throw std::invalid_argument("input data contains zeroes");
	}

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
					if (second_jet_it == med_btag_it)
					{
						continue;
					}

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

std::vector<int> to_one_hot(int pdgId)
{
	assert(abs(pdgId) == 13 || abs(pdgId) == 11);

	std::vector<int> one_hot;
	std::map<int, int> tmp{{-13, 0}, {-11, 0}, {11, 0}, {13, 0}};
	tmp[pdgId] += 1;

	std::transform(tmp.begin(), tmp.end(), std::back_inserter(one_hot), [](auto& pair){ return pair.second;});
	return one_hot;
}

int main() {
	// TFile *myFile = TFile::Open("out_radion_2L_2016_m400.root");
	TFile *myFile = TFile::Open("out_ttbar_2L_2016.root");
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

	float ak8_jet0_subjet1_E;
	myTree->SetBranchAddress("ak8_jet0_subjet1_E", &ak8_jet0_subjet1_E);

	float ak8_jet0_subjet1_px;
	myTree->SetBranchAddress("ak8_jet0_subjet1_px", &ak8_jet0_subjet1_px);

	float ak8_jet0_subjet1_py;
	myTree->SetBranchAddress("ak8_jet0_subjet1_py", &ak8_jet0_subjet1_py);

	float ak8_jet0_subjet1_pz;
	myTree->SetBranchAddress("ak8_jet0_subjet1_pz", &ak8_jet0_subjet1_pz);

	float ak8_jet0_subjet2_E;
	myTree->SetBranchAddress("ak8_jet0_subjet2_E", &ak8_jet0_subjet2_E);

	float ak8_jet0_subjet2_px;
	myTree->SetBranchAddress("ak8_jet0_subjet2_px", &ak8_jet0_subjet2_px);

	float ak8_jet0_subjet2_py;
	myTree->SetBranchAddress("ak8_jet0_subjet2_py", &ak8_jet0_subjet2_py);

	float ak8_jet0_subjet2_pz;
	myTree->SetBranchAddress("ak8_jet0_subjet2_pz", &ak8_jet0_subjet2_pz);

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

	float lep0_conept;
	myTree->SetBranchAddress("lep0_conept", &lep0_conept);

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

	float lep1_conept;
	myTree->SetBranchAddress("lep1_conept", &lep1_conept);

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

	float ak8_jet0_tau1;
	myTree->SetBranchAddress("ak8_jet0_tau1", &ak8_jet0_tau1);

	float ak8_jet0_tau2;
	myTree->SetBranchAddress("ak8_jet0_tau2", &ak8_jet0_tau2);

	float ak8_jet0_tau3;
	myTree->SetBranchAddress("ak8_jet0_tau3", &ak8_jet0_tau3);

	float ak8_jet0_tau4;
	myTree->SetBranchAddress("ak8_jet0_tau4", &ak8_jet0_tau4);

	float ak8_jet0_msoftdrop;
	myTree->SetBranchAddress("ak8_jet0_msoftdrop", &ak8_jet0_msoftdrop);

	int nEvents = myTree->GetEntries();
    auto branches = myTree->GetListOfBranches();
    int nBranches = myTree->GetNbranches();
	std::vector<std::string> branch_names;
	
	std::regex pat_jet("_jet[0-9]_[pEtm]");
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

	// std::copy(branch_names.begin(), branch_names.end(), std::ostream_iterator<std::string>(std::cout, "\n"));

	int same_bjets = 0;

	// the following 8 lines are needed to sort names of branches according to output of to_one_hot
	// map is ordered and sorted by key by default
	std::map<int, std::string> l0_map{{-13, "l0_anti_muon"}, {11, "l0_ele"}, {13, "l0_muon"}, {-11, "l0_anti_ele"}};
	std::map<int, std::string> l1_map{{-13, "l1_anti_muon"}, {11, "l1_ele"}, {13, "l1_muon"}, {-11, "l1_anti_ele"}};
	std::vector<std::string> one_hot_l0_names;
	std::vector<std::string> one_hot_l1_names;

	std::transform(l0_map.begin(), l0_map.end(), std::back_inserter(one_hot_l0_names), [](auto& pair){ return pair.second;});
	// std::copy(one_hot_l0_names.begin(), one_hot_l0_names.end(), std::ostream_iterator<std::string>(std::cout, ","));
	std::transform(l1_map.begin(), l1_map.end(), std::back_inserter(one_hot_l1_names), [](auto& pair){ return pair.second;});
	// std::copy(one_hot_l1_names.begin(), one_hot_l1_names.end(), std::ostream_iterator<std::string>(std::cout, ","));

	// std::vector<std::string> one_hot_l0_names{"l0_anti_muon", "l0_anti_ele", "l0_ele", "l0_muon"};
	// std::vector<std::string> one_hot_l1_names{"l1_anti_muon", "l1_anti_ele", "l1_ele", "l1_muon"};

	bool write_to_file = true;
	std::ofstream data;
	if (write_to_file)
	{
		// data.open("data_signal.csv");
		data.open("data_backgr.csv");

		// basic kinematics
		for (auto const& bName: branch_names)
		{
			if (bName == "lep0_pdgId" || bName == "lep1_pdgId")
			{
				continue;
			}
			data << bName << ",";
		} 

		// one-hot encoded leptons pdgId names
		std::ostream_iterator<std::string> out_it(data, ",");
		std::copy(one_hot_l0_names.begin(), one_hot_l0_names.end(), out_it);
		std::copy(one_hot_l1_names.begin(), one_hot_l1_names.end(), out_it);

		// high-level kinematics
		data << "ll_dR,"
			<< "ll_dPhi," 
			<< "ll_mass," 
			<< "llmet_dPhi,"
			<< "leadl_ak4jet_mindR,"
			<< "subl_ak4jet_mindR,"
			<< "min_jet_dR,"
			<< "min_jet_dPhi,"
			<< "lmet_mass,"
			<< "llmet_transverse_mass,"
			<< "llbbmet_mass," 
			<< "subj_l_mindR,"
			<< "leadj_l_mindR," 
			<< "llbb_dPhi,"
			<< "llbb_dR,"
			<< "bb_mass,"
			<< "bbmet_dPhi," 
			<< "bb_dPhi," 
			<< "bb_dR," 
			<< "lep0_conept," 
			<< "lep1_conept,"
			<< "n_btagged_ak4_jets,";
		data << "dnn_truth_value" << "\n";
	}

	std::map<int, int> btag_counts;
	btag_counts[0] = 0;
	btag_counts[1] = 0;
	btag_counts[2] = 0;
	btag_counts[3] = 0;
	btag_counts[4] = 0;

	std::map<int, int> zero_ak4jet_counts;
	zero_ak4jet_counts[0] = 0;
	zero_ak4jet_counts[1] = 0;
	zero_ak4jet_counts[2] = 0;
	zero_ak4jet_counts[3] = 0;
	zero_ak4jet_counts[4] = 0;

	std::map<int, int> nonzero_ak4jet_counts;
	nonzero_ak4jet_counts[0] = 0;
	nonzero_ak4jet_counts[1] = 0;
	nonzero_ak4jet_counts[2] = 0;
	nonzero_ak4jet_counts[3] = 0;
	nonzero_ak4jet_counts[4] = 0;

	std::map<int, int> zero_lep_couns;
	zero_lep_couns[0] = 0;
	zero_lep_couns[1] = 0;
	zero_lep_couns[2] = 0;

	std::map<std::string, int> jet_counts;
	jet_counts["boosted"] = 0;
	jet_counts["resolved"] = 0;

	int bad_jet_event_count = 0;
	int invalid_jet_event_count = 0;

	std::function<void(TLorentzVector const&)> print_4vector = [](TLorentzVector const& p) 
	{ 
		std::cout << "(" << p.E() << ", " 
				  << p.Px() << ", " 
				  << p.Py() << ", "
				  << p.Pz() << ")\n"; 
	};

	int n_skipped = 0;
	int n_equal_lep = 0;
	int n_equal_jet = 0;
	TLorentzVector const null(0.0, 0.0, 0.0, 0.0);
	std::cout << "nEvents = " << nEvents << std::endl;
	for (size_t i = 0; i < nEvents; ++i)
	{
		myTree->GetEntry(i);

		TLorentzVector l0, l1, met, ak4jet0, ak4jet1, ak4jet2, ak4jet3, ak8jet0, ak8jet0_subjet1, ak8jet0_subjet2;
		l0.SetPxPyPzE(lep0_px, lep0_py, lep0_pz, lep0_E);
		l1.SetPxPyPzE(lep1_px, lep1_py, lep1_pz, lep1_E);
		met.SetPxPyPzE(met_px, met_py, met_pz, met_E);
		ak4jet0.SetPxPyPzE(ak4_jet0_px, ak4_jet0_py, ak4_jet0_pz, ak4_jet0_E);
		ak4jet1.SetPxPyPzE(ak4_jet1_px, ak4_jet1_py, ak4_jet1_pz, ak4_jet1_E);
		ak4jet2.SetPxPyPzE(ak4_jet2_px, ak4_jet2_py, ak4_jet2_pz, ak4_jet2_E);
		ak4jet3.SetPxPyPzE(ak4_jet3_px, ak4_jet3_py, ak4_jet3_pz, ak4_jet3_E);

		ak8jet0.SetPxPyPzE(ak8_jet0_px, ak8_jet0_py, ak8_jet0_pz, ak8_jet0_E);
		ak8jet0_subjet1.SetPxPyPzE(ak8_jet0_subjet1_px, ak8_jet0_subjet1_py, ak8_jet0_subjet1_pz, ak8_jet0_subjet1_E);
		ak8jet0_subjet2.SetPxPyPzE(ak8_jet0_subjet2_px, ak8_jet0_subjet2_py, ak8_jet0_subjet2_pz, ak8_jet0_subjet2_E);

		std::vector<TLorentzVector> jets = {ak4jet0, ak4jet1, ak4jet2, ak4jet3};
		std::vector<TLorentzVector> leptons = {l0, l1};
		std::vector<float> btagScores = {ak4_jet0_btagDeepFlavB, ak4_jet1_btagDeepFlavB, ak4_jet2_btagDeepFlavB, ak4_jet3_btagDeepFlavB};
		std::multimap<float, TLorentzVector> jet_data;

		int n_zero_btags = std::count(btagScores.begin(), btagScores.end(), 0.0);
		btag_counts[4 - n_zero_btags] += 1;

		std::function<bool(TLorentzVector const&)> is_null = [&null](TLorentzVector const& p) { return p == null; };

		int n_zero_ak4jets = std::count_if(jets.begin(), jets.end(), is_null);
		zero_ak4jet_counts[n_zero_ak4jets] += 1;
		nonzero_ak4jet_counts[4 - n_zero_ak4jets] += 1;

		int n_zero_lep = std::count_if(leptons.begin(), leptons.end(), is_null);
		zero_lep_couns[n_zero_lep] += 1;

		bool resolved = (n_zero_ak4jets == 0 && n_zero_btags == 0);
		bool boosted = ((n_zero_ak4jets == 1 || n_zero_ak4jets == 2) && ak8jet0 != null && ak8jet0_subjet1 != null && ak8jet0_subjet2 != null);

		if (boosted)
		{
			jet_counts["boosted"] += 1;
			jets.push_back(ak8jet0_subjet1);
			jets.push_back(ak8jet0_subjet2);
			// std::copy(btagScores.begin(), btagScores.end(), std::ostream_iterator<float>(std::cout, " "));
			// std::cout << "\n";
		}
		

		if (resolved)
		{
			jet_counts["resolved"] += 1;
		}

		if (n_zero_lep == 0)
		{
			if (l0 == l1)
				++n_equal_lep;
		}

		// bool valid_event = (n_zero_jets < 3 && n_zero_btags < 3 && l0 != null && l1 != null);
		// bool valid_event = (n_zero_ak4jets == 0 && n_zero_btags == 0 && l0 != null && l1 != null && ak8jet0 != null);
		bool valid_event = (resolved || boosted) && l0 != null && l1 != null;

		if (!valid_event)
		{
			++n_skipped;
			continue;
		}

		for (size_t j = 0; j < jets.size(); ++j)
		{
			for (size_t k = j + 1; k < jets.size(); ++k)
			{
				if (jets[j] == jets[k] && jets[j] != null && jets[k] != null)
				{
					++n_equal_jet;
				}
			}
		}

		std::function<void(std::pair<float, TLorentzVector> const&)> print_jet_data = [](std::pair<float, TLorentzVector> const& p) 
		{ 
			std::cout << p.first << ": (" << p.second.E() << ", " 
					  					  << p.second.Px() << ", " 
										  << p.second.Py() << ", "
										  << p.second.Pz() << ")\n"; 
		};

		// insert ak4 jets into map with their btagging scores
		std::transform(btagScores.begin(), 
					   btagScores.end(), 
					   jets.begin(), 
					   std::inserter(jet_data, jet_data.end()), 
					   [](float a, TLorentzVector b) { return std::make_pair(a, b); });

		// in boosted event (2 or 3 ak4 jets and ak8 jet with 2 subjets) add subjets
		// the're assigned btagging score 1 because I don't know their btagging score and I want to preserve the original selection rule:
		// if there is a jet with medium btagging score, select the remaing jet so that invariant mass is closest to higgs mass
		// if there are no jets with medium btagging score, select the best pair by invariant mass criterion
		if (boosted)
		{
			jet_data.insert(std::make_pair(1.0, ak8jet0_subjet1));
			jet_data.insert(std::make_pair(1.0, ak8jet0_subjet2));
			// std::cout << "boosted:\n";
			// std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
			// std::cout << "-----------------------------------\n";
		}

		jet_data.erase(0); // now it contains only jets with nonzero momentum and btagging score

		std::vector<int> l0_one_hot = to_one_hot(lep0_pdgId);
		// std::cout << "event #" << i << ", lep0_pdgId = " << lep0_pdgId << "; l0_one_hot = ( ";
		// std::copy(l0_one_hot.begin(), l0_one_hot.end(), std::ostream_iterator<int>(std::cout, " "));
		// std::cout << ")" << std::endl;
		std::vector<int> l1_one_hot = to_one_hot(lep1_pdgId);
		// std::cout << "event #" << i << ", lep1_pdgId = " << lep1_pdgId << "; l1_one_hot = ( ";
		// std::copy(l1_one_hot.begin(), l1_one_hot.end(), std::ostream_iterator<int>(std::cout, " "));
		// std::cout << ")" << std::endl;
		// std::cout << "-------------------------------------------------\n";


		// for (auto const& jet: jet_data)
		// {
		// 	if (is_null(jet.second))
		// 	{
		// 		std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
		// 		std::cout << "-------------------------------------------------\n";
		// 	}
		// 	break;
		// }

		// if (std::find_if(jet_data.begin(), jet_data.end(), is_null))
		// {
		// 	std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
		// 	std::cout << "-------------------------------------------------\n";
		// }

		std::pair<TLorentzVector, TLorentzVector> jets_from_h;
		try 
		{
			jets_from_h = select_bjets(jet_data);
			float dijet_mass = (jets_from_h.first + jets_from_h.second).M();
			if (dijet_mass < 90.0 || dijet_mass > 160.0)
			{
				++bad_jet_event_count;
			}
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

		TLorentzVector lead_j, sub_j;
		if (valid_event)
		{
			if (jets_from_h.first.Pt() > jets_from_h.second.Pt())
			{
				lead_j = jets_from_h.first;
				sub_j = jets_from_h.second;
			}
			else
			{
				sub_j = jets_from_h.first;
				lead_j = jets_from_h.second;
			}
		}

		// if (is_null(lead_j) || is_null(sub_j))
		// {
		// 	std::cout << "jet from h is null in event " << i << std::endl;
		// }

		float bb_dR = lead_j.DeltaR(sub_j);
		float bb_dPhi = lead_j.DeltaPhi(sub_j);
		if (bb_dR == 0.0 || bb_dPhi == 0.0)
		{

			std::cout << "Event #" << i << ":" << std::endl;
			std::cout << "lead_j: ";
			print_4vector(lead_j);
			std::cout << "sub_j: ";
			print_4vector(sub_j);
			std::for_each(jet_data.begin(), jet_data.end(), print_jet_data);
			// std::for_each(jets.begin(), jets.end(), print_4vector);
			++same_bjets;
			std::cout << "---------------------------------------\n";
		}
		float bbmet_dPhi = met.DeltaPhi(lead_j + sub_j);
		float bb_mass = (lead_j + sub_j).M();

		float ll_dR = l0.DeltaR(l1);
		float ll_dPhi = l0.DeltaPhi(l1);
		float ll_mass = (l0 + l1).M();
		float llmet_dPhi = met.DeltaPhi(l0 + l1);
		
		TLorentzVector lead_l, sub_l;
		float lead_l_conept, sub_l_conept;
		if (l0.Pt() > l1.Pt())
		{
			lead_l = l0;
			sub_l = l1;
			lead_l_conept = lep0_conept;
			sub_l_conept = lep1_conept;
		}
		else
		{
			lead_l = l1;
			sub_l = l0;
			lead_l_conept = lep1_conept;
			sub_l_conept = lep0_conept;
		}

		float llbb_dR = (lead_j + sub_j).DeltaR(lead_l + sub_l);
		float llbb_dPhi = (lead_j + sub_j).DeltaPhi(lead_l + sub_l);

		auto dR = [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaR(p2); };
		auto dPhi = [](TLorentzVector const& p1, TLorentzVector const& p2) { return p1.DeltaPhi(p2); };

		// just in case, remove all zero vectors from jet collection
		jets.erase(std::remove_if(jets.begin(), jets.end(), is_null), jets.end());
		
		// if (std::count_if(jets.begin(), jets.end(), [](TLorentzVector const& p) { return p == TLorentzVector(0.0, 0.0, 0.0, 0.0); }) != 0)
		// {
		// 	std::cout << "Event " << i << " contains zero vectors:\n";
		// 	std::for_each(jets.begin(), jets.end(), print_4vector);
		// 	std::cout << "---------------------------------------\n";
		// }

		// for (auto const& jet: jets)
		// {
		// 	if (is_null(jet))
		// 	{
		// 		std::for_each(jets.begin(), jets.end(), print_4vector);
		// 		std::cout << "---------------------------------------\n";
		// 		break;
		// 	}
		// }

		// if (boosted)
		// {
		// 	std::for_each(jets.begin(), jets.end(), print_4vector);
		// 	std::cout << "---------------------------------------\n";
		// }

		float leadl_ak4jet_mindR = min_func(lead_l, jets, dR);
		float subl_ak4jet_mindR = min_func(sub_l, jets, dR);

		float leadj_l_mindR = min_func(lead_j, leptons, dR);
		float subj_l_mindR = min_func(sub_j, leptons, dR);

		float min_jet_dR = min_func(jets, dR);
		float min_jet_dPhi = min_func(jets, dPhi);

		float lmet_mass = (lead_l + met).M();
		float llbbmet_mass = (lead_l + sub_l + sub_j + lead_j + met).M();

		float llmet_transverse_mass;
		float dphi_llmet = TVector2::Phi_mpi_pi((l0 + l1).Phi() - met.Phi());
		llmet_transverse_mass = sqrt(2.0*(l0 + l1).Pt()*met.Pt()*(1 - cos(dphi_llmet)));

		int n_btagged_ak4_jets;
		if (resolved)
		{
			n_btagged_ak4_jets = jet_data.size();
		}
		if (boosted)
		{
			n_btagged_ak4_jets = jet_data.size() - 2;
		}
		if (n_btagged_ak4_jets <= 0)
		{
			std::cout << "Event #" << i << ": n_btagged_ak4_jets <= 0\n";
		}

		// std::cout << "llmet_transverse_mass = " << llmet_transverse_mass << std::endl;

		// std::cout << "Event: " << i << std::endl
		// 		  << "ll_dR = " << ll_dR << "\n" 
		// 		  << "ll_dPhi = " << ll_dPhi << "\n"
		// 		  << "ll_mass = " << ll_mass << "\n"
		// 		  << "llmet_dPhi = " << llmet_dPhi << "\n"
		// 		  << "leadl_ak4jet_mindR = " << leadl_ak4jet_mindR << "\n"
		// 		  << "subl_ak4jet_mindR = " << subl_ak4jet_mindR << "\n"
		// 		  << "min_jet_dR = " << min_jet_dR << "\n"
		// 		  << "min_jet_dPhi = " << min_jet_dPhi << "\n"
		// 		  << "lmet_mass = " << lmet_mass << "\n"
		// 		  << "llmet_transverse_mass = " << llmet_transverse_mass << "\n"
		// 		  << "llbbmet_mass = " << llbbmet_mass << "\n"
		// 		  << "subj_l_mindR = " << subj_l_mindR << "\n"
		// 		  << "leadj_l_mindR = " << leadj_l_mindR << "\n"
		// 		  << "llbb_dPhi = " << llbb_dPhi << "\n"
		// 		  << "llbb_dR = " << llbb_dR << "\n"
		// 		  << "bb_mass = " << bb_mass << "\n"
		// 		  << "bbmet_dPhi = " << bbmet_dPhi << "\n"
		// 		  << "bb_dPhi = " << bb_dPhi << "\n"
		// 		  << "bb_dR = " << bb_dR << "\n"
		// 		  << "lead_l_conept = " << lead_l_conept << "\n"
		// 		  << "sub_l_conept = " << sub_l_conept << "\n"
		// 		  << "--------------------------------\n";
	
		// break;

		if (write_to_file)
		{
			// basic kinematics
			data << ak4_jet0_E << ","
				 << ak4_jet0_phi << ","
				 << ak4_jet0_pt << ","
				 << ak4_jet0_px << ","
				 << ak4_jet0_py << ","
				 << ak4_jet0_pz << ","
				 << ak4_jet1_E << ","
				 << ak4_jet1_phi << ","
				 << ak4_jet1_pt << ","
				 << ak4_jet1_px << ","
				 << ak4_jet1_py << ","
				 << ak4_jet1_pz << ","
				 << ak4_jet2_E << ","
				 << ak4_jet2_phi << ","
				 << ak4_jet2_pt << ","
				 << ak4_jet2_px << ","
				 << ak4_jet2_py << ","
				 << ak4_jet2_pz << ","
				 << ak4_jet3_E << ","
				 << ak4_jet3_phi << ","
				 << ak4_jet3_pt << ","
				 << ak4_jet3_px << ","
				 << ak4_jet3_py << ","
				 << ak4_jet3_pz << ","
				 << ak8_jet0_E << ","
				 << ak8_jet0_msoftdrop << ","
				 << ak8_jet0_phi << ","
				 << ak8_jet0_pt << ","
				 << ak8_jet0_px << ","
				 << ak8_jet0_py << ","
				 << ak8_jet0_pz << ","
				 << ak8_jet0_tau1 << ","
				 << ak8_jet0_tau2 << ","
				 << ak8_jet0_tau3 << ","
				 << ak8_jet0_tau4 << ","
				 << lep0_E << ","
				 << lep0_phi << ","
				 << lep0_pt << ","
				 << lep0_px << ","
				 << lep0_py << ","
				 << lep0_pz << ","
				 << lep1_E << ","
				 << lep1_phi << ","
				 << lep1_pt << ","
				 << lep1_px << ","
				 << lep1_py << ","
				 << lep1_pz << ","
				 << met_E << ","
				 << met_px << ","
				 << met_py << ","
				 << met_pz << ",";

			// one-hot encoded leptons pdgId
			std::copy(l0_one_hot.begin(), l0_one_hot.end(), std::ostream_iterator<int>(data, ","));
			std::copy(l1_one_hot.begin(), l1_one_hot.end(), std::ostream_iterator<int>(data, ","));

			// high-level kinematics
			data << ll_dR << "," 
				 << ll_dPhi << ","
				 << ll_mass << ","
				 << llmet_dPhi << ","
				 << leadl_ak4jet_mindR << ","
				 << subl_ak4jet_mindR << ","
				 << min_jet_dR << ","
				 << min_jet_dPhi << ","
				 << lmet_mass << ","
				 << llmet_transverse_mass << ","
				 << llbbmet_mass << ","
				 << subj_l_mindR << ","
				 << leadj_l_mindR << ","
				 << llbb_dPhi << ","
				 << llbb_dR << ","
				 << bb_mass << ","
				 << bbmet_dPhi << ","
			 	 << bb_dPhi << ","
				 << bb_dR << ","
				 << lep0_conept << ","
				 << lep1_conept << ","
				 << n_btagged_ak4_jets << ","
				 << dnn_truth_value << "\n";
		}
		// break;
	}

	std::cout << "---------------------------\n";
	std::for_each(btag_counts.begin(), btag_counts.end(), [](std::pair<int, int> const& p) { std::cout << p.first << " b-tagged jet(s): " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::for_each(zero_ak4jet_counts.begin(), zero_ak4jet_counts.end(), [](std::pair<int, int> const& p) { std::cout << p.first << " zero jet(s): " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::for_each(zero_lep_couns.begin(), zero_lep_couns.end(), [](std::pair<int, int> const& p) { std::cout << p.first << " zero lepton(s): " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::for_each(jet_counts.begin(), jet_counts.end(), [](std::pair<std::string, int> const& p) { std::cout << p.first << ": " << p.second << std::endl; });
	std::cout << "---------------------------\n";
	std::cout << "b jets not selected: invalid_jet_event_count = " << invalid_jet_event_count << std::endl;
	std::cout << "higgs from bb mass wrong: bad_jet_event_count = " << bad_jet_event_count << std::endl;
	std::cout << "Total events skipped n_skipped = " << n_skipped << std::endl;
	std::cout << "Total events with same b jets = " << same_bjets << std::endl;
	std::cout << "Total events with equal leptons n_equal_lep = " << n_equal_lep << std::endl;
	std::cout << "Total events with equal jets n_equal_jet = " << n_equal_jet << std::endl;
	std::cout << "---------------------------\n";
	
	if (write_to_file)
	{
		data.close();
	}

	return 0;
}
