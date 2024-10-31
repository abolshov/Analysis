#include "Validator.hpp"
#include "EstimatorTools.hpp"

#include <iostream>

#include "TCanvas.h"
#include "TLine.h"

Validator::Validator()
{   
    rg.SetSeed(42);

    // momentum ratios of quarks and jets
    hm.Add("corr_ratio_px_1", "Ratio of corrected Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_py_1", "Ratio of corrected Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_pz_1", "Ratio of corrected Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_E_1", "Ratio of corrected E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("corr_ratio_px_2", "Ratio of corrected Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_py_2", "Ratio of corrected Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_pz_2", "Ratio of corrected Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_E_2", "Ratio of corrected E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("raw_ratio_px_1", "Ratio of Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_py_1", "Ratio of Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_pz_1", "Ratio of Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_E_1", "Ratio of E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("raw_ratio_px_2", "Ratio of Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_py_2", "Ratio of Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_pz_2", "Ratio of Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_E_2", "Ratio of E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 25}, 100);

    // momenta of quarks and jets
    hm.Add("corr_px_1", "Corrected Px of jet 1", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_py_1", "Corrected Py of jet 1", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_pz_1", "Corrected Pz of jet 1", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_E_1", "Corrected E of jet 1", {"E, [GeV]", "Count"}, {0.0, 400.0}, 100);

    hm.Add("corr_px_2", "Corrected Px of jet 2", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_py_2", "Corrected Py of jet 2", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_pz_2", "Corrected Pz of jet 2", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_E_2", "Corrected E of jet 2", {"E, [GeV]", "Count"}, {0.0, 400.0}, 100);

    hm.Add("raw_px_1", "Raw Px of jet 1", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_py_1", "Raw Py of jet 1", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_pz_1", "Raw Pz of jet 1", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_E_1", "Raw E of jet 1", {"E, [GeV]", "Count"}, {0.0, 400.0}, 100);

    hm.Add("raw_px_2", "Raw Px of jet 2", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_py_2", "Raw Py of jet 2", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_pz_2", "Raw Pz of jet 2", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_E_2", "Raw E of jet 2", {"E, [GeV]", "Count"}, {0.0, 400.0}, 100);

    hm.Add("true_px_1", "True Px of quark 1", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_py_1", "True Py of quark 1", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_pz_1", "True Pz of quark 1", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_E_1", "True E of quark 1", {"E, [GeV]", "Count"}, {0, 400.0}, 100);

    hm.Add("true_px_2", "True Px of quark 2", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_py_2", "True Py of quark 2", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_pz_2", "True Pz of quark 2", {"Pz, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_E_2", "True E of quark 2", {"E, [GeV]", "Count"}, {0, 400.0}, 100);

    // px and py of met
    hm.Add("true_px_met", "True (gen) Px of MET", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_py_met", "True (gen) Py of MET", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);

    hm.Add("corr_px_met", "Corrected Px of MET", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("corr_py_met", "Corrected Py of MET", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);

    hm.Add("raw_px_met", "Raw Px of MET", {"Px, [GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("raw_py_met", "Raw Py of MET", {"Py, [GeV]", "Count"}, {-400.0, 400.0}, 100);

    // differences
    hm.Add("true_corr_px_1_diff", "Difference between corrected Px and true Px of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_py_1_diff", "Difference between corrected Py and true Py of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_pz_1_diff", "Difference between corrected Pz and true Pz of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_E_1_diff", "Difference between corrected E and true E of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);

    hm.Add("true_corr_px_2_diff", "Difference between corrected Px and true Px of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_py_2_diff", "Difference between corrected Py and true Py of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_pz_2_diff", "Difference between corrected Pz and true Pz of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_corr_E_2_diff", "Difference between corrected E and true E of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);

    hm.Add("true_raw_px_1_diff", "Difference between raw Px and true Px of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_py_1_diff", "Difference between raw Py and true Py of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_pz_1_diff", "Difference between raw Pz and true Pz of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_E_1_diff", "Difference between raw E and true E of jet 1", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);

    hm.Add("true_raw_px_2_diff", "Difference between raw Px and true Px of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_py_2_diff", "Difference between raw Py and true Py of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_pz_2_diff", "Difference between raw Pz and true Pz of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
    hm.Add("true_raw_E_2_diff", "Difference between raw E and true E of jet 2", {"[GeV]", "Count"}, {-400.0, 400.0}, 100);
}

void Validator::FillVariables(Event const& event)
{
    recoMET = LorentzVectorF_t(event.reco_met_pt, 0.0, event.reco_met_phi, 0.0);

    for (int i = 0; i < event.reco_jet.nRecoJet; ++i)
    {
        reco_jet_p4.emplace_back(event.reco_jet.pt[i], event.reco_jet.eta[i], event.reco_jet.phi[i], event.reco_jet.mass[i]);

        jet_PNet_resolutions.push_back(event.reco_jet.pt[i]*event.reco_jet.PNetRegPtRawRes[i]);
        jet_PNet_corrections.push_back(event.reco_jet.PNetRegPtRawCorr[i]);
    }

    // gen truth inputs
    gen_truth_p4.resize(event.m_index.size(), LorentzVectorF_t{});
    for (auto const& [obj_name, idx]: event.m_index)
    {
        gen_truth_p4[idx] = LorentzVectorF_t(event.gen_truth.pt[idx], event.gen_truth.eta[idx], event.gen_truth.phi[idx], event.gen_truth.mass[idx]);
    }

    // true nu inputs
    nu.resize(event.m_nu_index.size(), LorentzVectorF_t{});
    for (auto const& [obj_name, idx]: event.m_nu_index)
    {
        nu[idx] = LorentzVectorF_t(event.nu.pt[idx], event.nu.eta[idx], event.nu.phi[idx], event.nu.mass[idx]);
    }
}

void Validator::ClearVariables()
{
    gen_truth_p4.clear();
    reco_jet_p4.clear(); 
    nu.clear();
    jet_PNet_resolutions.clear(); 
    jet_PNet_corrections.clear();
}

void Validator::Compare(long long event_number, bool use_2d_pdf, bool skip_failures)
{
    // reorder b1 and b2 to make quark with larger pt b1
    if (gen_truth_p4[truth_index["b1"]].Pt() < gen_truth_p4[truth_index["b2"]].Pt())
    {
        std::swap(gen_truth_p4[truth_index["b1"]], gen_truth_p4[truth_index["b2"]]);
    }

    LorentzVectorF_t const& bq1 = gen_truth_p4[truth_index["b1"]];
    LorentzVectorF_t const& bq2 = gen_truth_p4[truth_index["b2"]];

    LorentzVectorF_t raw_bj1, raw_bj2;
    if (reco_jet_p4[0].Pt() > reco_jet_p4[1].Pt())
    {
        raw_bj1 = reco_jet_p4[0];
        raw_bj2 = reco_jet_p4[1];
    }
    else
    {
        raw_bj2 = reco_jet_p4[0];
        raw_bj1 = reco_jet_p4[1];
    }

    LorentzVectorF_t const& genMET = gen_truth_p4[truth_index["met"]];

    hm.Fill("true_px_met", genMET.Px());
    hm.Fill("true_py_met", genMET.Py());

    double reco_met_px = recoMET.Px();
    double reco_met_py = recoMET.Py();

    hm.Fill("true_px_1", bq1.Px());
    hm.Fill("true_py_1", bq1.Py());
    hm.Fill("true_pz_1", bq1.Pz());
    hm.Fill("true_E_1", bq1.E());

    hm.Fill("true_px_2", bq2.Px());
    hm.Fill("true_py_2", bq2.Py());
    hm.Fill("true_pz_2", bq2.Pz());
    hm.Fill("true_E_2", bq2.E());

    UHist1d_t& pdf1d = pdf_1d[pdf1d_index["pdf_b1"]];
    UHist1d_t& pdf_mbb = pdf_1d[pdf1d_index["pdf_mbb"]];
    UHist2d_t& pdf2d = pdf_2d[pdf2d_index["pdf_b1b2"]];

    // histograms for accumulating effective ratio values for one event
    auto hist_ratio_px_1 = std::make_unique<TH1F>("hist_ratio_px_1", "hist_ratio_px_1", 100, -25, 25);
    auto hist_ratio_py_1 = std::make_unique<TH1F>("hist_ratio_py_1", "hist_ratio_py_1", 100, -25, 25);
    auto hist_ratio_pz_1 = std::make_unique<TH1F>("hist_ratio_pz_1", "hist_ratio_pz_1", 100, -25, 25);
    auto hist_ratio_E_1 = std::make_unique<TH1F>("hist_ratio_E_1", "hist_ratio_E_1", 100, 0, 25);

    auto hist_ratio_px_2 = std::make_unique<TH1F>("hist_ratio_px_2", "hist_ratio_px_2", 100, -25, 25);
    auto hist_ratio_py_2 = std::make_unique<TH1F>("hist_ratio_py_2", "hist_ratio_py_2", 100, -25, 25);
    auto hist_ratio_pz_2 = std::make_unique<TH1F>("hist_ratio_pz_2", "hist_ratio_pz_2", 100, -25, 25);
    auto hist_ratio_E_2 = std::make_unique<TH1F>("hist_ratio_E_2", "hist_ratio_E_2", 100, 0, 25);

    // histograms for accumulating effective momentum values for one event
    auto hist_px_1 = std::make_unique<TH1F>("hist_px_1", "hist_px_1", 100, -400.0, 400.0);
    auto hist_py_1 = std::make_unique<TH1F>("hist_py_1", "hist_py_1", 100, -400.0, 400.0);
    auto hist_pz_1 = std::make_unique<TH1F>("hist_pz_1", "hist_pz_1", 100, -400.0, 400.0);
    auto hist_E_1 = std::make_unique<TH1F>("hist_E_1", "hist_E_1", 100, 0.0, 400.0);

    auto hist_px_2 = std::make_unique<TH1F>("hist_px_2", "hist_px_2", 100, -400.0, 400.0);
    auto hist_py_2 = std::make_unique<TH1F>("hist_py_2", "hist_py_2", 100, -400.0, 400.0);
    auto hist_pz_2 = std::make_unique<TH1F>("hist_pz_2", "hist_pz_2", 100, -400.0, 400.0);
    auto hist_E_2 = std::make_unique<TH1F>("hist_E_2", "hist_E_2", 100, 0.0, 400.0);

    // histograms for accumulating effective MET momentum values for one event
    auto hist_px_met = std::make_unique<TH1F>("hist_px_met", "hist_px_met", 100, -400.0, 400.0);
    auto hist_py_met = std::make_unique<TH1F>("hist_py_met", "hist_py_met", 100, -400.0, 400.0);

    hm.Fill("raw_ratio_px_1", raw_bj1.Px()/bq1.Px());
    hm.Fill("raw_ratio_py_1", raw_bj1.Py()/bq1.Py());
    hm.Fill("raw_ratio_pz_1", raw_bj1.Pz()/bq1.Pz());
    hm.Fill("raw_ratio_E_1", raw_bj1.E()/bq1.E());

    hm.Fill("raw_ratio_px_2", raw_bj2.Px()/bq2.Px());
    hm.Fill("raw_ratio_py_2", raw_bj2.Py()/bq2.Py());
    hm.Fill("raw_ratio_pz_2", raw_bj2.Pz()/bq2.Pz());
    hm.Fill("raw_ratio_E_2", raw_bj2.E()/bq2.E());

    hm.Fill("raw_px_1", raw_bj1.Px());
    hm.Fill("raw_py_1", raw_bj1.Py());
    hm.Fill("raw_pz_1", raw_bj1.Pz());
    hm.Fill("raw_E_1", raw_bj1.E());

    hm.Fill("raw_px_2", raw_bj2.Px());
    hm.Fill("raw_py_2", raw_bj2.Py());
    hm.Fill("raw_pz_2", raw_bj2.Pz());
    hm.Fill("raw_E_2", raw_bj2.E());

    hm.Fill("raw_px_met", reco_met_px);
    hm.Fill("raw_py_met", reco_met_py);

    for (int i = 0; i < N_ITER; ++i)
    {
        LorentzVectorF_t bj1 = raw_bj1;
        LorentzVectorF_t bj2 = raw_bj2;

        double c1 = 1.0;
        double c2 = 1.0;

        if (use_2d_pdf)
        {
            pdf2d->GetRandom2(c1, c2, &rg);
        }
        else 
        {
            auto resc_pair = JetRescFact(bj1, bj2, pdf1d, HIGGS_MASS);
            if (resc_pair)
            {
                auto [resc1, resc2] = resc_pair.value();
                c1 = resc1;
                c2 = resc2;
            }
            else 
            {
                if (skip_failures)
                {
                    continue;
                }    
            }
        }

        bj1 *= c1;
        bj2 *= c2;

        double w = 1.0;
        if (use_2d_pdf)
        {
            double mbb = (bj1 + bj2).M();
            int bin = pdf_mbb->FindBin(mbb);
            w = pdf_mbb->GetBinContent(bin);
        }

        // corrections to met due to b jets
        double dpx_bjet = -(c1 - 1)*bj1.Px() - (c2 - 1)*bj2.Px();
        double dpy_bjet = -(c1 - 1)*bj1.Py() - (c2 - 1)*bj2.Py();

        // std::cout << "dpx_bjet=" << dpx_bjet << ", dpy_bjet=" << dpy_bjet << "\n";

        // corrections to met due to smearing
        double dpx_smear = rg.Gaus(0, MET_SIGMA);
        double dpy_smear = rg.Gaus(0, MET_SIGMA);
        
        double reco_met_corr_px = reco_met_px + dpx_bjet + dpx_smear;
        double reco_met_corr_py = reco_met_py + dpy_bjet + dpy_smear;

        if (reco_met_corr_px > 400.0)
        {
            std::cout << "event " << event_number << ":\n";
            std::cout << "\tc1=" << c1 << ", c2=" << c2 << "\n";

            std::cout << "\traw_bj1=(" << raw_bj1.E() << ", " << raw_bj1.Px() << ", " << raw_bj1.Py() << ", " << raw_bj1.Pz() << ")\n";
            std::cout << "\tcor_bj1=(" << bj1.E() << ", " << bj1.Px() << ", " << bj1.Py() << ", " << bj1.Pz() << ")\n";  

            std::cout << "\traw_bj2=(" << raw_bj2.E() << ", " << raw_bj2.Px() << ", " << raw_bj2.Py() << ", " << raw_bj2.Pz() << ")\n";
            std::cout << "\tcor_bj2=(" << bj2.E() << ", " << bj2.Px() << ", " << bj2.Py() << ", " << bj2.Pz() << ")\n";  

            std::cout << "\tm(bb)=" << (bj1 + bj2).M() << "\n";
            std::cout << "\tw=" << w << "\n";

            std::cout << "\treco_met_corr_px=" << reco_met_corr_px 
                      << ", reco_met_px=" << reco_met_px 
                      << ", dpx_bjet=" << dpx_bjet
                      << ", dpx_smear=" << dpx_smear << "\n";

            std::cout << "\treco_met_corr_py=" << reco_met_corr_py 
                      << ", reco_met_py=" << reco_met_py 
                      << ", dpy_bjet=" << dpy_bjet
                      << ", dpy_smear=" << dpy_smear << "\n";
            std::cout << "\n";

            std::cin.get();
        }

        hist_px_met->Fill(reco_met_corr_px);
        hist_py_met->Fill(reco_met_corr_py);

        hist_px_1->Fill(bj1.Px(), w);
        hist_py_1->Fill(bj1.Py(), w);
        hist_pz_1->Fill(bj1.Pz(), w);
        hist_E_1->Fill(bj1.E(), w);

        hist_px_2->Fill(bj2.Px(), w);
        hist_py_2->Fill(bj2.Py(), w);
        hist_pz_2->Fill(bj2.Pz(), w);
        hist_E_2->Fill(bj2.E(), w);

        double ratio_px_1 = bj1.Px()/bq1.Px();
        double ratio_py_1 = bj1.Py()/bq1.Py();
        double ratio_pz_1 = bj1.Pz()/bq1.Pz();
        double ratio_E_1 = bj1.E()/bq1.E();

        hist_ratio_px_1->Fill(ratio_px_1, w);
        hist_ratio_py_1->Fill(ratio_py_1, w);
        hist_ratio_pz_1->Fill(ratio_pz_1, w);
        hist_ratio_E_1->Fill(ratio_E_1, w);

        double ratio_px_2 = bj2.Px()/bq2.Px();
        double ratio_py_2 = bj2.Py()/bq2.Py();
        double ratio_pz_2 = bj2.Pz()/bq2.Pz();
        double ratio_E_2 = bj2.E()/bq2.E();

        hist_ratio_px_2->Fill(ratio_px_2, w);
        hist_ratio_py_2->Fill(ratio_py_2, w);
        hist_ratio_pz_2->Fill(ratio_pz_2, w);
        hist_ratio_E_2->Fill(ratio_E_2, w);
    }

    auto GetEffValue = [](UHist1d_t& h){ return h->GetEntries() ? h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) : 0.0; };

    if (event_number % 413 == 54)
    {
        auto c1 = std::make_unique<TCanvas>("c1", "c1");

        char const* pdf_type = use_2d_pdf ? "pdf2d" : "pdf1d";

        hist_px_1->Draw();
        double x1 = bq1.Px();
        double y2 = hist_px_1->GetBinContent(hist_px_1->GetMaximumBin());
        TLine l_true_px_1(x1, 0.0, x1, y2);
        l_true_px_1.SetLineColor(kRed);
        l_true_px_1.Draw("same");
        TLine l_raw_px_1(raw_bj1.Px(), 0.0, raw_bj1.Px(), y2);
        l_raw_px_1.SetLineColor(kGreen);
        l_raw_px_1.Draw("same");
        c1->SaveAs(Form("valid_plots/px/px1_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_px_2->Draw();
        x1 = bq2.Px();
        y2 = hist_px_2->GetBinContent(hist_px_2->GetMaximumBin());
        TLine l_true_px_2(x1, 0.0, x1, y2);
        l_true_px_2.SetLineColor(kRed);
        l_true_px_2.Draw("same");
        TLine l_raw_px_2(raw_bj2.Px(), 0.0, raw_bj2.Px(), y2);
        l_raw_px_2.SetLineColor(kGreen);
        l_raw_px_2.Draw("same");
        c1->SaveAs(Form("valid_plots/px/px2_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_py_1->Draw();
        x1 = bq1.Py();
        y2 = hist_py_1->GetBinContent(hist_py_1->GetMaximumBin());
        TLine l_true_py_1(x1, 0.0, x1, y2);
        l_true_py_1.SetLineColor(kRed);
        l_true_py_1.Draw("same");
        TLine l_raw_py_1(raw_bj1.Py(), 0.0, raw_bj1.Py(), y2);
        l_raw_py_1.SetLineColor(kGreen);
        l_raw_py_1.Draw("same");
        c1->SaveAs(Form("valid_plots/py/py1_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_py_2->Draw();
        x1 = bq2.Py();
        y2 = hist_py_2->GetBinContent(hist_py_2->GetMaximumBin());
        TLine l_true_py_2(x1, 0.0, x1, y2);
        l_true_py_2.SetLineColor(kRed);
        l_true_py_2.Draw("same");
        TLine l_raw_py_2(raw_bj2.Py(), 0.0, raw_bj2.Py(), y2);
        l_raw_py_2.SetLineColor(kGreen);
        l_raw_py_2.Draw("same");
        c1->SaveAs(Form("valid_plots/py/py2_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_pz_1->Draw();
        x1 = bq1.Pz();
        y2 = hist_pz_1->GetBinContent(hist_pz_1->GetMaximumBin());
        TLine l_true_pz_1(x1, 0.0, x1, y2);
        l_true_pz_1.SetLineColor(kRed);
        l_true_pz_1.Draw("same");
        TLine l_raw_pz_1(raw_bj1.Pz(), 0.0, raw_bj1.Pz(), y2);
        l_raw_pz_1.SetLineColor(kGreen);
        l_raw_pz_1.Draw("same");
        c1->SaveAs(Form("valid_plots/pz/pz1_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_pz_2->Draw();
        x1 = bq2.Pz();
        y2 = hist_pz_2->GetBinContent(hist_pz_2->GetMaximumBin());
        TLine l_true_pz_2(x1, 0.0, x1, y2);
        l_true_pz_2.SetLineColor(kRed);
        l_true_pz_2.Draw("same");
        TLine l_raw_pz_2(raw_bj2.Pz(), 0.0, raw_bj2.Pz(), y2);
        l_raw_pz_2.SetLineColor(kGreen);
        l_raw_pz_2.Draw("same");
        c1->SaveAs(Form("valid_plots/pz/pz2_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_E_1->Draw();
        x1 = bq1.E();
        y2 = hist_E_1->GetBinContent(hist_E_1->GetMaximumBin());
        TLine l_true_E_1(x1, 0.0, x1, y2);
        l_true_E_1.SetLineColor(kRed);
        l_true_E_1.Draw("same");
        TLine l_raw_E_1(raw_bj1.E(), 0.0, raw_bj1.E(), y2);
        l_raw_E_1.SetLineColor(kGreen);
        l_raw_E_1.Draw("same");
        c1->SaveAs(Form("valid_plots/E/E1_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));

        hist_E_2->Draw();
        x1 = bq2.E();
        y2 = hist_E_2->GetBinContent(hist_E_2->GetMaximumBin());
        TLine l_true_E_2(x1, 0.0, x1, y2);
        l_true_E_2.SetLineColor(kRed);
        l_true_E_2.Draw("same");
        TLine l_raw_E_2(raw_bj2.E(), 0.0, raw_bj2.E(), y2);
        l_raw_E_2.SetLineColor(kGreen);
        l_raw_E_2.Draw("same");
        c1->SaveAs(Form("valid_plots/E/E2_%s_evt_%d.png", pdf_type, static_cast<int>(event_number)));
    }
    
    double eff_ratio_px_1 = GetEffValue(hist_ratio_px_1);
    double eff_ratio_py_1 = GetEffValue(hist_ratio_py_1);
    double eff_ratio_pz_1 = GetEffValue(hist_ratio_pz_1);
    double eff_ratio_E_1 = GetEffValue(hist_ratio_E_1);

    double eff_ratio_px_2 = GetEffValue(hist_ratio_px_2);
    double eff_ratio_py_2 = GetEffValue(hist_ratio_py_2);
    double eff_ratio_pz_2 = GetEffValue(hist_ratio_pz_2);
    double eff_ratio_E_2 = GetEffValue(hist_ratio_E_2);

    hm.Fill("corr_ratio_px_1", eff_ratio_px_1);
    hm.Fill("corr_ratio_py_1", eff_ratio_py_1);
    hm.Fill("corr_ratio_pz_1", eff_ratio_pz_1);
    hm.Fill("corr_ratio_E_1", eff_ratio_E_1);

    hm.Fill("corr_ratio_px_2", eff_ratio_px_2);
    hm.Fill("corr_ratio_py_2", eff_ratio_py_2);
    hm.Fill("corr_ratio_pz_2", eff_ratio_pz_2);
    hm.Fill("corr_ratio_E_2", eff_ratio_E_2);

    double eff_px_1 = GetEffValue(hist_px_1);
    double eff_py_1 = GetEffValue(hist_py_1);
    double eff_pz_1 = GetEffValue(hist_pz_1);
    double eff_E_1 = GetEffValue(hist_E_1);

    double eff_px_2 = GetEffValue(hist_px_2);
    double eff_py_2 = GetEffValue(hist_py_2);
    double eff_pz_2 = GetEffValue(hist_pz_2);
    double eff_E_2 = GetEffValue(hist_E_2);

    hm.Fill("corr_px_1", eff_px_1);
    hm.Fill("corr_py_1", eff_py_1);
    hm.Fill("corr_pz_1", eff_pz_1);
    hm.Fill("corr_E_1", eff_E_1);

    hm.Fill("corr_px_2", eff_px_2);
    hm.Fill("corr_py_2", eff_py_2);
    hm.Fill("corr_pz_2", eff_pz_2);
    hm.Fill("corr_E_2", eff_E_2);

    double eff_px_met = GetEffValue(hist_px_met);
    double eff_py_met = GetEffValue(hist_py_met);

    // std::cout << "reco_met_px=" << reco_met_px << ", eff_px_met=" << eff_px_met << ", true_met_px=" << genMET.Px() << "\n";

    hm.Fill("corr_px_met", eff_px_met);
    hm.Fill("corr_py_met", eff_py_met);

    hm.Fill("true_corr_px_1_diff", bq1.Px() - eff_px_1);
    hm.Fill("true_corr_py_1_diff", bq1.Py() - eff_py_1);
    hm.Fill("true_corr_pz_1_diff", bq1.Pz() - eff_pz_1);
    hm.Fill("true_corr_E_1_diff", bq1.E() - eff_E_1);

    hm.Fill("true_corr_px_2_diff", bq2.Px() - eff_px_2);
    hm.Fill("true_corr_py_2_diff", bq2.Py() - eff_py_2);
    hm.Fill("true_corr_pz_2_diff", bq2.Pz() - eff_pz_2);
    hm.Fill("true_corr_E_2_diff", bq2.E() - eff_E_2);

    hm.Fill("true_raw_px_1_diff", bq1.Px() - raw_bj1.Px());
    hm.Fill("true_raw_py_1_diff", bq1.Py() - raw_bj1.Py());
    hm.Fill("true_raw_pz_1_diff", bq1.Pz() - raw_bj1.Pz());
    hm.Fill("true_raw_E_1_diff", bq1.E() - raw_bj1.E());

    hm.Fill("true_raw_px_2_diff", bq2.Px() - raw_bj2.Px());
    hm.Fill("true_raw_py_2_diff", bq2.Py() - raw_bj2.Py());
    hm.Fill("true_raw_pz_2_diff", bq2.Pz() - raw_bj2.Pz());
    hm.Fill("true_raw_E_2_diff", bq2.E() - raw_bj2.E());

    ClearVariables();
}
