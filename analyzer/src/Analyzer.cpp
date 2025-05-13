#include "Analyzer.hpp"
#include "Definitions.hpp"
#include "EstimatorUtils.hpp"

#ifdef DEV
    #include "MatchingTools.hpp"
    #include "SelectionUtils.hpp"
#endif 

#include <iostream>
#include <chrono>

#include "TH1.h"

Analyzer::Analyzer(std::map<Channel, TString> const& pdf_file_map)
:   m_estimator(pdf_file_map.at(Channel::SL), pdf_file_map.at(Channel::DL), AggregationMode::Event)
,   m_estimator_data(std::make_unique<EstimatorData>())
{
    TH1::AddDirectory(false);
}

void Analyzer::ProcessDataset(Dataset_t const& dataset)
{
    for (auto const& sample: dataset)
    {
        ProcessSample(sample);
    }
}

void Analyzer::ProcessSample(Sample const& sample)
{
    counter = 0;

    TString const& file_name = sample.file_name;
    TString const& tree_name = sample.tree_name;
    TString const& type = sample.type;
    Channel ch = sample.channel;
    Int_t masspoint = sample.masspoint;
    bool is_bkg = masspoint < 0;

    std::cout << "Start processing " << file_name << "\n";
    auto process_file_start{std::chrono::steady_clock::now()};
    auto input_file = std::make_unique<TFile>(file_name);
    auto input_tree = std::unique_ptr<TTree>(static_cast<TTree*>(input_file->Get<TTree>(tree_name)));

    if (m_record_output)
    {
        TString output_file_name = FormFileName(ch, is_bkg, masspoint, type, "out");
        std::cout << "Output data will be recorded to " << output_file_name << "\n";
        m_recorder.OpenFile(output_file_name);
        m_recorder.ResetTree(MakeTree());
    }

    UTree_t& output_tree = m_recorder.GetTree(); // be careful, may be ref to nullptr

    #ifdef DEV
        if (m_record_iterations)
        {
            TString iter_file_name = FormFileName(ch, is_bkg, masspoint, type, "iter");
            std::cout << "Iteration data will be recorded to " << iter_file_name << "\n";
            m_estimator.OpenDbgFile(iter_file_name, ch);
        }
    #endif

    m_event.ConnectTree(input_tree, ch);
    ULong64_t n_events = input_tree->GetEntries();
    std::chrono::duration<double> average_duration{};
    for (ULong64_t evt = 0; evt < n_events; ++evt)
    {
        auto process_event_start{std::chrono::steady_clock::now()};
        #ifdef DEV 
            if (m_process_matched_slice && !is_bkg)
            {
                ProcessEventSlice(evt, input_tree, output_tree, ch);
            }
            else 
            {
                ProcessEvent(evt, input_tree, output_tree, ch, is_bkg);
            }
        #else
            ProcessEvent(evt, input_tree, output_tree, ch, is_bkg);
        #endif
        auto process_event_end{std::chrono::steady_clock::now()};
        std::chrono::duration<double> process_event_duration{process_event_end - process_event_start};
        average_duration += process_event_duration;
        if (evt % 5000 == 0 && evt > 0)
        {
            std::cout << "===> Processed " << 100.0*evt/n_events << "%, time elapsed " 
                      << std::chrono::duration<double>(average_duration).count() << " s, average time per event: " 
                      << std::chrono::duration<double, std::milli>(average_duration).count()/counter << " ms\n";
        }
    }
    std::cout << "counter=" << counter << "\n";

    if (m_record_output)
    {
        m_recorder.WriteTree(OUT_TREE_NAME);
    }

    // is needed to avoid double free 
    // if done immediately after getting tree, root delets tree cache and nothing is written to tree
    input_tree->SetDirectory(nullptr); 
    input_file->Close();

    auto process_file_end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds{process_file_end - process_file_start};
    std::cout << "Processing " << file_name << " took " << elapsed_seconds.count() << " seconds\n";
    std::cout << "Average time per event: " << std::chrono::duration<double, std::milli>(average_duration).count()/counter << " ms\n";
}

void Analyzer::ProcessEvent(ULong64_t evt, UTree_t& input_tree, UTree_t& output_tree, Channel ch, bool is_bkg)
{
    input_tree->GetEntry(evt);

    ULong64_t evt_id = m_event.event_id;
    if (evt_id % 2 != 1)
    {
        return;
    }

    VecLVF_t jets = GetRecoJetP4(m_event);
    // VecLVF_t leptons = GetRecoLepP4(m_event, ch);
    // LorentzVectorF_t met = GetRecoMET(m_event);
    
    #ifdef DEV 
        if (!is_bkg)
        {
            if (!IsRecoverable(m_event, ch))
            {
                return;
            }

            if (!IsFiducial(m_event, jets, ch))
            {
                return;
            }
        }
    #endif

    ++counter;
    m_estimator_data->event_id = evt_id;

    // std::vector<Float_t> resolutions = GetPNetRes(m_event);
    // auto hme = m_estimator.EstimateMass(jets, resolutions, leptons, met, evt_id, ch);
    auto hme = m_estimator.EstimateMass(m_event, ch);
    if (hme)
    {
        auto out_array = hme.value();
        m_estimator_data->mass = out_array[static_cast<size_t>(EstimOut::mass)];
        m_estimator_data->integral = out_array[static_cast<size_t>(EstimOut::integral)];
        m_estimator_data->peak_value = out_array[static_cast<size_t>(EstimOut::peak_value)];
        m_estimator_data->width = out_array[static_cast<size_t>(EstimOut::width)];
    }
    else 
    {
        m_estimator_data->mass = -1.0f;
        m_estimator_data->integral = -1.0f;
        m_estimator_data->peak_value = -1.0f;
        m_estimator_data->width = -1.0f;
    }

    if (m_record_output && output_tree != nullptr)
    {
        m_recorder.FillTree();
    }
}

UTree_t Analyzer::MakeTree()
{
    auto output_tree = std::make_unique<TTree>(OUT_TREE_NAME, "HME output tree");
    output_tree->SetDirectory(nullptr); 

    output_tree->Branch("hme_mass", &m_estimator_data->mass, "hme_mass/F");
    output_tree->Branch("hme_integral", &m_estimator_data->integral, "hme_integral/F");
    output_tree->Branch("hme_width", &m_estimator_data->width, "hme_width/F");
    output_tree->Branch("hme_peak_value", &m_estimator_data->peak_value, "hme_peak_value/F");
    output_tree->Branch("event_id", &m_estimator_data->event_id, "event_id/l");

    return output_tree;
}

#ifdef DEV 
    void Analyzer::ProcessEventSlice(ULong64_t evt, UTree_t& input_tree, UTree_t& output_tree, Channel ch)
    {
        input_tree->GetEntry(evt);

        VecLVF_t jets = GetRecoJetP4(m_event);
        VecLVF_t quarks = GetGenQuarksP4(m_event, ch);
        VecLVF_t leptons = GetRecoLepP4(m_event, ch);
        LorentzVectorF_t met = GetRecoMET(m_event);

        ULong64_t event_id = m_event.event_id;
        if (event_id % 2 != 1)
        {
            return;
        }

        if (!IsRecoverable(m_event, ch))
        {
            return;
        }

        if (!IsFiducial(m_event, jets, ch))
        {
            return;
        }

        JetComb match = FindMatch(quarks, jets, ch);
        if (!match.HasUniqueJets(ch))
        {
            return;
        }

        ++counter;
        m_estimator_data->event_id = event_id;
        std::vector<Float_t> resolutions = GetPNetRes(m_event);
        
        VecLVF_t particles;
        if (ch == Channel::SL)
        {
            particles.resize(static_cast<size_t>(ObjSL::count));
            particles[static_cast<size_t>(ObjSL::lep)] = leptons[static_cast<size_t>(Lep::lep1)];
            particles[static_cast<size_t>(ObjSL::met)] = met;    
            
            if (jets[match.b1].Pt() > jets[match.b2].Pt())
            {
                particles[static_cast<size_t>(ObjSL::bj1)] = jets[match.b1];
                particles[static_cast<size_t>(ObjSL::bj2)] = jets[match.b2];
            }
            else 
            {
                particles[static_cast<size_t>(ObjSL::bj2)] = jets[match.b1];
                particles[static_cast<size_t>(ObjSL::bj1)] = jets[match.b2];
            }

            if (jets[match.q1].Pt() > jets[match.q2].Pt())
            {
                particles[static_cast<size_t>(ObjSL::lj1)] = jets[match.q1];
                particles[static_cast<size_t>(ObjSL::lj2)] = jets[match.q2];
            }
            else 
            {
                particles[static_cast<size_t>(ObjSL::lj2)] = jets[match.q1];
                particles[static_cast<size_t>(ObjSL::lj1)] = jets[match.q2];
            }
        }
        else if (ch == Channel::DL)
        {
            particles.resize(static_cast<size_t>(ObjDL::count));
            particles[static_cast<size_t>(ObjDL::lep1)] = leptons[static_cast<size_t>(Lep::lep1)];
            particles[static_cast<size_t>(ObjDL::lep2)] = leptons[static_cast<size_t>(Lep::lep2)];
            particles[static_cast<size_t>(ObjDL::met)] = met;    
            if (jets[match.b1].Pt() > jets[match.b2].Pt())
            {
                particles[static_cast<size_t>(ObjDL::bj1)] = jets[match.b1];
                particles[static_cast<size_t>(ObjDL::bj2)] = jets[match.b2];
            }
            else 
            {
                particles[static_cast<size_t>(ObjDL::bj2)] = jets[match.b1];
                particles[static_cast<size_t>(ObjDL::bj1)] = jets[match.b2];
            }
        }
        else 
        {
            throw std::runtime_error("Unknown channel");
        }

        ArrF_t<ESTIM_OUT_SZ> comb_result = m_estimator.GetEstimator(ch).EstimateCombination(particles, resolutions, event_id, match);
        if (comb_result[static_cast<size_t>(EstimOut::mass)] > 0.0)
        {
            m_estimator_data->mass = comb_result[static_cast<size_t>(EstimOut::mass)];
            m_estimator_data->integral = comb_result[static_cast<size_t>(EstimOut::integral)];
            m_estimator_data->peak_value = comb_result[static_cast<size_t>(EstimOut::peak_value)];
            m_estimator_data->width = comb_result[static_cast<size_t>(EstimOut::width)];
        }
        else 
        {
            m_estimator_data->mass = -1.0f;
            m_estimator_data->integral = -1.0f;
            m_estimator_data->peak_value = -1.0f;
            m_estimator_data->width = -1.0f;
        }

        if (m_record_output && output_tree != nullptr)
        {
            m_recorder.FillTree();
        }
    }
#endif

TString Analyzer::FormFileName(Channel ch, bool is_bkg, Int_t mp, TString const& sample_type, TString const& file_type) const
{
    TString file_name;
    if (!is_bkg) 
    {
        #ifdef DEV
            if (m_process_matched_slice)
            {
                file_name = ch == Channel::SL ? Form("hme_%s_true_sl_M%d.root", file_type.Data(), mp) : Form("hme_%s_true_dl_M%d.root", file_type.Data(), mp);
            }
            else 
            {
                file_name = ch == Channel::SL ? Form("hme_%s_sl_M%d.root", file_type.Data(), mp) : Form("hme_%s_dl_M%d.root", file_type.Data(), mp);
            }
        #else 
            file_name = ch == Channel::SL ? Form("hme_%s_sl_M%d.root", file_type.Data(), mp) : Form("hme_%s_dl_M%d.root", file_type.Data(), mp);
        #endif
    }
    else
    {
        file_name = ch == Channel::SL ? Form("hme_%s_sl_%s.root", file_type.Data(), sample_type.Data()) : Form("hme_%s_dl_%s.root", file_type.Data(), sample_type.Data());
    }
    return file_name;
}