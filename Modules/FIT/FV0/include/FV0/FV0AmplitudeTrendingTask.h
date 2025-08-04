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
/// \file   FV0AmplitudeTrendingTask.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Custom trending task for FV0 amplitude analysis with auto-discovery
///

#ifndef QC_MODULE_FV0_AMPLITUDETRENDINGTASK_H
#define QC_MODULE_FV0_AMPLITUDETRENDINGTASK_H

#include "QualityControl/PostProcessingInterface.h"
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <array>

class TAxis;
class TH1;

namespace o2::quality_control::repository
{
class DatabaseInterface;
}

namespace o2::quality_control_modules::fv0
{

/// FV0 detector ring definition
struct FV0DetectorRing {
    std::string name;           ///< Ring name (R1, R2, R3, R4, R51, R52)
    std::string displayName;    ///< Display name (A1-H1, A2-H2, etc.)
    int startChannel;           ///< First channel in ring
    int endChannel;             ///< Last channel in ring
    int color;                  ///< Ring-specific color
    
    FV0DetectorRing() = default;
    FV0DetectorRing(const std::string& n, const std::string& dn, int start, int end, int col)
        : name(n), displayName(dn), startChannel(start), endChannel(end), color(col) {}
    
    std::vector<int> getChannels() const {
        std::vector<int> channels;
        for (int ch = startChannel; ch <= endChannel; ++ch) {
            channels.push_back(ch);
        }
        return channels;
    }
    
    std::string getChannelName(int channel) const {
        if (channel < startChannel || channel > endChannel) return "";
        int sectorIndex = channel - startChannel;
        char sectorLetter = 'A' + sectorIndex;
        return std::string(1, sectorLetter) + name.substr(1); // A1, B1, etc.
    }
};

/// Configuration for beam type
struct BeamConfig {
    std::string name;           ///< Short name (pp, PbPb, OO)
    std::string displayName;    ///< Display name (pp collisions, etc.)
    double expectedGain;        ///< Expected gain value for scaling
    
    BeamConfig() = default;
    BeamConfig(const std::string& n, const std::string& dn, double gain = 15.0)
        : name(n), displayName(dn), expectedGain(gain) {}
};

/// Plot template configuration
struct PlotTemplate {
    std::string name;           ///< Template name
    std::string type;           ///< Type: "per_channel", "per_ring", "cross_config", "all_channels"
    std::string title;          ///< Plot title template
    bool enabled;               ///< Whether to generate this template
    
    PlotTemplate() = default;
    PlotTemplate(const std::string& n, const std::string& t, const std::string& title_template, bool en = true)
        : name(n), type(t), title(title_template), enabled(en) {}
};

/// Custom trending task for FV0 amplitude analysis
class FV0AmplitudeTrendingTask final : public quality_control::postprocessing::PostProcessingInterface
{
public:
    FV0AmplitudeTrendingTask() = default;
    ~FV0AmplitudeTrendingTask() override = default;
    
    void configure(const boost::property_tree::ptree& config) override;
    void initialize(quality_control::postprocessing::Trigger trigger,
                   framework::ServiceRegistryRef services) override;
    void update(quality_control::postprocessing::Trigger trigger,
               framework::ServiceRegistryRef services) override;
    void finalize(quality_control::postprocessing::Trigger trigger,
                 framework::ServiceRegistryRef services) override;

private:
    // Configuration
    std::string mScalarsPath{"TrendsScalars"};          ///< Base path for trending scalars
    std::string mSourceTaskName{"AmplitudePostProc"};   ///< Source task name
    std::vector<BeamConfig> mBeamConfigs;               ///< Available beam configurations
    std::vector<PlotTemplate> mPlotTemplates;           ///< Plot templates to generate
    static constexpr int sNCHANNELS = 48;               ///< Number of FV0 channels
    
    // FV0 detector structure
    std::vector<FV0DetectorRing> mDetectorRings;
    
    // Trending data
    std::unique_ptr<TTree> mTrend;
    std::map<std::string, std::unique_ptr<TObject>> mPlots;
    
    // Metadata for trending
    struct {
        Long64_t runNumber = 0;
        static const char* getBranchLeafList() { return "runNumber/L"; }
    } mMetaData;
    UInt_t mTime;
    
    // Channel data arrays (one per beam config)
    std::map<std::string, std::vector<double>> mChannelData;
    
    // Configuration and initialization
    void initializeDetectorStructure();
    void configureBeamConfigs(const boost::property_tree::ptree& config);
    void configurePlotTemplates(const boost::property_tree::ptree& config);
    void autoDiscoverAvailableConfigs(quality_control::repository::DatabaseInterface& qcdb);
    void initializeTrending(quality_control::repository::DatabaseInterface& qcdb);
    
    // Data collection
    bool collectTrendingData(const quality_control::postprocessing::Trigger& t, 
                            quality_control::repository::DatabaseInterface& qcdb);
    
    // Plot generation
    void generateAllPlots();
    void generatePerChannelPlots();
    void generatePerRingPlots();  
    void generateCrossConfigPlots();
    void generateAllChannelsPlots();
    
    // Utility methods
    TCanvas* createBasicCanvas(const std::string& name, const std::string& title);
    void applyChannelStyling(TGraphErrors* graph, int channel, const std::string& beamConfig);
    void setupAxes(TH1* hist, const std::string& title, const std::string& xLabel, const std::string& yLabel);
    void formatTimeAxis(TH1* hist);
    std::string getChannelDisplayName(int channel) const;
    const FV0DetectorRing* getRingForChannel(int channel) const;
    
    // Constants for styling
    static constexpr std::array<int, 6> sRingColors = {44, 46, 47, 38, 29, 28};     // One color per ring
    static constexpr std::array<int, 8> sMarkerStyles = {21, 22, 23, 24, 25, 26, 27, 28}; // One per sector
};

} // namespace o2::quality_control_modules::fv0

#endif // QC_MODULE_FV0_AMPLITUDETRENDINGTASK_H