// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   FV0AmplitudeTrendingTask.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Custom trending task for FV0 amplitude analysis with auto-discovery
///

#include "FV0/FV0AmplitudeTrendingTask.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/DatabaseInterface.h"
#include "QualityControl/RepoPathUtils.h"

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace o2::quality_control;
using namespace o2::quality_control::core;
using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control::repository;

namespace o2::quality_control_modules::fv0
{

void FV0AmplitudeTrendingTask::configure(const boost::property_tree::ptree& config)
{
    const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
    const char* cfgCustom = Form("%s.custom", cfgPath);
    
    auto cfgGet = [&cfgCustom](const std::string& entry) {
        return Form("%s.%s", cfgCustom, entry.c_str());
    };
    
    // Basic configuration
    mScalarsPath = config.get<std::string>(cfgGet("scalarsPath"), "TrendsScalars");
    mSourceTaskName = config.get<std::string>(cfgGet("sourceTaskName"), "AmplitudePostProc");
    
    // Configure beam configurations and plot templates
    configureBeamConfigs(config);
    configurePlotTemplates(config);
    initializeDetectorStructure();
    
    ILOG(Info, Support) << "FV0AmplitudeTrendingTask configured with " << mBeamConfigs.size() 
                        << " beam configs and " << mPlotTemplates.size() << " plot templates" << ENDM;
}

void FV0AmplitudeTrendingTask::configureBeamConfigs(const boost::property_tree::ptree& config)
{
    const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
    const char* cfgCustom = Form("%s.custom", cfgPath);
    
    mBeamConfigs.clear();
    
    try {
        auto beamConfigNode = config.get_child(Form("%s.beamConfigs", cfgCustom));
        for (const auto& configItem : beamConfigNode) {
            std::string name = configItem.second.get<std::string>("name");
            std::string displayName = configItem.second.get<std::string>("displayName", name);
            double expectedGain = configItem.second.get<double>("expectedGain", 15.0);
            
            mBeamConfigs.emplace_back(name, displayName, expectedGain);
            ILOG(Info, Support) << "Configured beam config: " << name << " (" << displayName 
                               << ") with expected gain " << expectedGain << ENDM;
        }
    } catch (const std::exception& e) {
        // Default configurations
        mBeamConfigs.emplace_back("pp", "pp collisions", 15.0);
        mBeamConfigs.emplace_back("PbPb", "Pb-Pb collisions", 4.0);
        mBeamConfigs.emplace_back("OO", "O-O collisions", 11.0);
        ILOG(Info, Support) << "Using default beam configurations" << ENDM;
    }
}

void FV0AmplitudeTrendingTask::configurePlotTemplates(const boost::property_tree::ptree& config)
{
    const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
    const char* cfgCustom = Form("%s.custom", cfgPath);
    
    mPlotTemplates.clear();
    
    try {
        auto templateNode = config.get_child(Form("%s.plotTemplates", cfgCustom));
        for (const auto& templateItem : templateNode) {
            std::string name = templateItem.second.get<std::string>("name");
            std::string type = templateItem.second.get<std::string>("type");
            std::string title = templateItem.second.get<std::string>("title", name);
            bool enabled = templateItem.second.get<bool>("enabled", true);
            
            mPlotTemplates.emplace_back(name, type, title, enabled);
        }
    } catch (const std::exception& e) {
        // Default plot templates
        mPlotTemplates.emplace_back("per_ring", "per_ring", "FV0 {ring} - {config} Gain Trends", true);
        mPlotTemplates.emplace_back("cross_config", "cross_config", "FV0 {channel} - Cross Configuration", true);
        mPlotTemplates.emplace_back("inner_vs_outer", "ring_comparison", "FV0 Inner vs Outer Rings", true);
        ILOG(Info, Support) << "Using default plot templates" << ENDM;
    }
}

void FV0AmplitudeTrendingTask::initializeDetectorStructure()
{
    mDetectorRings.clear();
    
    // Define FV0 ring structure with colors
    mDetectorRings.emplace_back("R1", "A1-H1", 0, 7, sRingColors[0]);
    mDetectorRings.emplace_back("R2", "A2-H2", 8, 15, sRingColors[1]);
    mDetectorRings.emplace_back("R3", "A3-H3", 16, 23, sRingColors[2]);
    mDetectorRings.emplace_back("R4", "A4-H4", 24, 31, sRingColors[3]);
    mDetectorRings.emplace_back("R51", "A51-H51", 32, 39, sRingColors[4]);
    mDetectorRings.emplace_back("R52", "A52-H52", 40, 47, sRingColors[5]);
    
    ILOG(Info, Support) << "Initialized FV0 detector structure with " << mDetectorRings.size() << " rings" << ENDM;
}

void FV0AmplitudeTrendingTask::initializeTrending(DatabaseInterface& qcdb)
{
    // Create new TTree for trending
    mTrend = std::make_unique<TTree>();
    mTrend->SetName(PostProcessingInterface::getName().c_str());
    
    // Add meta and time branches
    mTrend->Branch("meta", &mMetaData, mMetaData.getBranchLeafList());
    mTrend->Branch("time", &mTime);
    
    // Add branches for each available beam config and channel
    for (const auto& beamConfig : mBeamConfigs) {
        mChannelData[beamConfig.name].resize(sNCHANNELS, 0.0);
        
        for (int ch = 0; ch < sNCHANNELS; ++ch) {
            std::string branchName = beamConfig.name + Form("_ch%02d", ch);
            mTrend->Branch(branchName.c_str(), &mChannelData[beamConfig.name][ch]);
        }
    }
    
    ILOG(Info, Support) << "Initialized trending TTree with " << mTrend->GetNbranches() << " branches" << ENDM;
}

void FV0AmplitudeTrendingTask::initialize(Trigger trigger, framework::ServiceRegistryRef services)
{
    auto& qcdb = services.get<DatabaseInterface>();
    
    // Initialize trending
    initializeTrending(qcdb);
    
    // Clear any existing plots
    mPlots.clear();
    getObjectsManager()->startPublishing(mTrend.get(), PublicationPolicy::ThroughStop);
    
    ILOG(Info, Support) << "FV0AmplitudeTrendingTask initialized successfully" << ENDM;
}

bool FV0AmplitudeTrendingTask::collectTrendingData(const Trigger& t, DatabaseInterface& qcdb)
{
  // Use producer's validFrom for x-axis; read with "latest valid for t.timestamp"
  const unsigned long long validFromMs =
      static_cast<unsigned long long>(t.activity.mValidity.getMin() ? t.activity.mValidity.getMin()
                                                                    : t.timestamp);
  mTime = validFromMs / 1000; // ROOT expects seconds
  mMetaData.runNumber = t.activity.mId;

  bool seenAny = false;

  for (const auto& cfg : mBeamConfigs) {
    for (int ch = 0; ch < sNCHANNELS; ++ch) {
      // Base folder + object name
      const std::string folder =
          RepoPathUtils::getMoPath("FV0",                      /*detector*/
                                   mSourceTaskName.c_str(),    /*task*/
                                   (mScalarsPath + "/" + cfg.name).c_str(),
                                   /*subsub=*/"", /*withQcPrefix*/ false);
      const std::string objName = cfg.name + Form("_ch%02d", ch);

      // 1) Primary: latest object valid for the trigger timestamp
      auto mo = qcdb.retrieveMO(folder, objName, DatabaseInterface::Timestamp::Latest);

      // 2) Fallback: explicit validFrom subfolder
      if (!mo || !mo->getObject()) {
        const std::string folderWithTs = folder + "/" + std::to_string(validFromMs);
        mo = qcdb.retrieveMO(folderWithTs, objName, DatabaseInterface::Timestamp::Latest);
      }

      if (mo && mo->getObject()) {
        if (auto* h = dynamic_cast<TH1F*>(mo->getObject()); h && h->GetEntries() > 0) {
          mChannelData[cfg.name][ch] = h->GetMean();
          seenAny = true;
        } else {
          mChannelData[cfg.name][ch] = 0.0;
        }
      } else {
        // Leave zero; just means this channel wasn’t there at this time
        mChannelData[cfg.name][ch] = 0.0;
      }
    }
  }

  return seenAny; // Fill only if at least one channel existed
}

void FV0AmplitudeTrendingTask::update(Trigger trigger, framework::ServiceRegistryRef services)
{
    auto& qcdb = services.get<DatabaseInterface>();
    
    // Collect trending data
    bool dataCollected = collectTrendingData(trigger, qcdb);
    
    if (dataCollected) {
        mTrend->Fill();
        generateAllPlots();
    } else {
        ILOG(Warning, Support) << "Not all data was available for trending" << ENDM;
    }
}

void FV0AmplitudeTrendingTask::finalize(Trigger, framework::ServiceRegistryRef)
{
    // Publish the trending tree
    getObjectsManager()->startPublishing(mTrend.get());
    
    // Generate final plots
    generateAllPlots();
    
    ILOG(Info, Support) << "FV0AmplitudeTrendingTask finalized with " << mTrend->GetEntries() 
                        << " trend entries and " << mPlots.size() << " plots" << ENDM;
}

void FV0AmplitudeTrendingTask::generateAllPlots()
{
    if (!mTrend || mTrend->GetEntries() < 1) {
        ILOG(Info, Support) << "No trending data available, skipping plot generation" << ENDM;
        return;
    }
    
    ILOG(Info, Support) << "Generating plots from " << mTrend->GetEntries() << " trend entries..." << ENDM;
    
    for (const auto& plotTemplate : mPlotTemplates) {
        if (!plotTemplate.enabled) continue;
        
        if (plotTemplate.type == "per_ring") {
            generatePerRingPlots();
        } else if (plotTemplate.type == "cross_config") {
            generateCrossConfigPlots();
        } else if (plotTemplate.type == "ring_comparison") {
            generateAllChannelsPlots();
        }
    }
}

void FV0AmplitudeTrendingTask::generatePerRingPlots()
{
    for (const auto& ring : mDetectorRings) {
        for (const auto& beamConfig : mBeamConfigs) {
            
            std::string plotName = Form("PerRing/%s_%s_ring_trend", beamConfig.name.c_str(), ring.name.c_str());
            std::string plotTitle = Form("FV0 %s - %s Gain Trends", ring.displayName.c_str(), beamConfig.displayName.c_str());
            
            auto* canvas = createBasicCanvas(plotName, plotTitle);
            auto* legend = new TLegend(0.75, 0.65, 0.95, 0.95);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            
            bool firstGraph = true;
            auto channels = ring.getChannels();
            
            for (size_t i = 0; i < channels.size(); ++i) {
                int ch = channels[i];
                std::string varexp = Form("%s_ch%02d:time", beamConfig.name.c_str(), ch);
                std::string drawOption = firstGraph ? "AL*" : "SAME L*";
                
                mTrend->Draw(varexp.c_str(), "", drawOption.c_str());
                
                if (auto* graph = dynamic_cast<TGraphErrors*>(canvas->FindObject("Graph"))) {
                    applyChannelStyling(graph, ch, beamConfig.name);
                    graph->SetName(Form("ch%02d", ch));
                    
                    std::string channelName = ring.getChannelName(ch);
                    legend->AddEntry(graph, channelName.c_str(), "lp");
                }
                
                firstGraph = false;
            }
            
            // Set up axes
            if (auto* htemp = dynamic_cast<TH1*>(canvas->FindObject("htemp"))) {
                setupAxes(htemp, plotTitle, "Time", "Gain (ADC/MIP)");
                formatTimeAxis(htemp);
                htemp->GetYaxis()->SetRangeUser(0.5 * beamConfig.expectedGain, 1.5 * beamConfig.expectedGain);
            }
            
            legend->Draw();
            canvas->Modified();
            canvas->Update();
            
            mPlots[plotName] = std::unique_ptr<TObject>(canvas);
            getObjectsManager()->startPublishing(canvas, PublicationPolicy::Once);
        }
    }
}

void FV0AmplitudeTrendingTask::generateCrossConfigPlots()
{
    // Generate plots for selected channels across all beam configs
    std::vector<int> selectedChannels = {0, 8, 16, 24, 32, 40}; // One from each ring
    
    for (int ch : selectedChannels) {
        std::string plotName = Form("CrossConfig/ch%02d_cross_config_trend", ch);
        std::string plotTitle = Form("FV0 %s - Cross Configuration Comparison", getChannelDisplayName(ch).c_str());
        
        auto* canvas = createBasicCanvas(plotName, plotTitle);
        auto* legend = new TLegend(0.75, 0.75, 0.95, 0.95);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        
        bool firstGraph = true;
        
        for (const auto& beamConfig : mBeamConfigs) {
            
            std::string varexp = Form("%s_ch%02d:time", beamConfig.name.c_str(), ch);
            std::string drawOption = firstGraph ? "AL*" : "SAME L*";
            
            mTrend->Draw(varexp.c_str(), "", drawOption.c_str());
            
            if (auto* graph = dynamic_cast<TGraphErrors*>(canvas->FindObject("Graph"))) {
                applyChannelStyling(graph, ch, beamConfig.name);
                graph->SetName(beamConfig.name.c_str());
                legend->AddEntry(graph, beamConfig.displayName.c_str(), "lp");
            }
            
            firstGraph = false;
        }
        
        // Set up axes
        if (auto* htemp = dynamic_cast<TH1*>(canvas->FindObject("htemp"))) {
            setupAxes(htemp, plotTitle, "Time", "Gain (ADC/MIP)");
            formatTimeAxis(htemp);
            htemp->GetYaxis()->SetRangeUser(2.0, 18.0); // Wide range for cross-config
        }
        
        legend->Draw();
        canvas->Modified();
        canvas->Update();
        
        mPlots[plotName] = std::unique_ptr<TObject>(canvas);
        getObjectsManager()->startPublishing(canvas, PublicationPolicy::Once);
    }
}

void FV0AmplitudeTrendingTask::generateAllChannelsPlots()
{
    // Generate summary plots showing ring averages
    for (const auto& beamConfig : mBeamConfigs) {
        
        std::string plotName = Form("Summary/%s_ring_summary", beamConfig.name.c_str());
        std::string plotTitle = Form("FV0 %s - Ring Summary Trends", beamConfig.displayName.c_str());
        
        auto* canvas = createBasicCanvas(plotName, plotTitle);
        auto* legend = new TLegend(0.75, 0.70, 0.95, 0.95);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        
        // For now, just plot representative channels from each ring
        bool firstGraph = true;
        for (const auto& ring : mDetectorRings) {
            int representativeChannel = ring.startChannel; // First channel of each ring
            
            std::string varexp = Form("%s_ch%02d:time", beamConfig.name.c_str(), representativeChannel);
            std::string drawOption = firstGraph ? "AL*" : "SAME L*";
            
            mTrend->Draw(varexp.c_str(), "", drawOption.c_str());
            
            if (auto* graph = dynamic_cast<TGraphErrors*>(canvas->FindObject("Graph"))) {
                graph->SetLineColor(ring.color);
                graph->SetMarkerColor(ring.color);
                graph->SetMarkerStyle(21);
                graph->SetMarkerSize(1.2);
                graph->SetLineWidth(2);
                graph->SetName(ring.name.c_str());
                
                legend->AddEntry(graph, ring.displayName.c_str(), "lp");
            }
            
            firstGraph = false;
        }
        
        // Set up axes
        if (auto* htemp = dynamic_cast<TH1*>(canvas->FindObject("htemp"))) {
            setupAxes(htemp, plotTitle, "Time", "Gain (ADC/MIP)");
            formatTimeAxis(htemp);
            htemp->GetYaxis()->SetRangeUser(0.5 * beamConfig.expectedGain, 1.5 * beamConfig.expectedGain);
        }
        
        legend->Draw();
        canvas->Modified();
        canvas->Update();
        
        mPlots[plotName] = std::unique_ptr<TObject>(canvas);
        getObjectsManager()->startPublishing(canvas, PublicationPolicy::Once);
    }
}

TCanvas* FV0AmplitudeTrendingTask::createBasicCanvas(const std::string& name, const std::string& title)
{
    auto* canvas = new TCanvas(name.c_str(), title.c_str(), 1200, 800);
    canvas->SetGrid();
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.12);
    return canvas;
}

void FV0AmplitudeTrendingTask::applyChannelStyling(TGraphErrors* graph, int channel, const std::string& beamConfig)
{
    if (!graph) return;
    
    const FV0DetectorRing* ring = getRingForChannel(channel);
    if (!ring) return;
    
    int sectorIndex = channel - ring->startChannel;
    
    graph->SetLineColor(ring->color);
    graph->SetMarkerColor(ring->color);
    graph->SetMarkerStyle(sMarkerStyles[sectorIndex % sMarkerStyles.size()]);
    graph->SetMarkerSize(1.1);
    graph->SetLineWidth(2);
}

void FV0AmplitudeTrendingTask::setupAxes(TH1* hist, const std::string& title, 
                                        const std::string& xLabel, const std::string& yLabel)
{
    if (!hist) return;
    
    hist->SetTitle(title.c_str());
    hist->GetXaxis()->SetTitle(xLabel.c_str());
    hist->GetYaxis()->SetTitle(yLabel.c_str());
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
}

void FV0AmplitudeTrendingTask::formatTimeAxis(TH1* hist)
{
    if (!hist) return;
    
    hist->GetXaxis()->SetTimeDisplay(1);
    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetTimeOffset(0.0);
    hist->GetXaxis()->SetTimeFormat("%m-%d %H:%M");
}

std::string FV0AmplitudeTrendingTask::getChannelDisplayName(int channel) const
{
    const FV0DetectorRing* ring = getRingForChannel(channel);
    if (!ring) return Form("Ch%02d", channel);
    
    return ring->getChannelName(channel);
}

const FV0DetectorRing* FV0AmplitudeTrendingTask::getRingForChannel(int channel) const
{
    for (const auto& ring : mDetectorRings) {
        if (channel >= ring.startChannel && channel <= ring.endChannel) {
            return &ring;
        }
    }
    return nullptr;
}

} // namespace o2::quality_control_modules::fv0