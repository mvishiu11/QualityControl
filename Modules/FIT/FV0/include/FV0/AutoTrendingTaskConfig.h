#ifndef QUALITYCONTROL_AUTOTRENDINGTASKCONFIG_H
#define QUALITYCONTROL_AUTOTRENDINGTASKCONFIG_H

#include "QualityControl/TrendingTaskConfig.h"
#include <boost/property_tree/ptree_fwd.hpp>
#include <string>
#include <vector>

namespace o2::quality_control_modules::fv0
{

/// Lightweight knobs to auto-generate FV0 trending configuration.
struct AutoTrendingConfig {
  // Trending behavior
  bool        producePlotsOnUpdate{true};
  bool        resumeTrend{true};
  bool        trendIfAllInputs{false};
  std::string trendingTimestamp{"validFrom"};

  /// Beam configuration (e.g. pp, OO, PbPb).
  struct BeamConfig {
    std::string name;        ///< key used in data paths (e.g. "pp")
    std::string displayName; ///< human-friendly label
    double      expectedGain{15.0}; ///< nominal gain g; y-range is [0.5g, 1.5g]
  };

  /// Channel layout and data source.
  struct ChannelConfig {
    int         totalChannels{48};
    int         channelsPerRing{8};
    std::string dataSourcePath{"FV0/MO/AmplitudePostProc/TrendsScalars"};
    std::string reductorName{"o2::quality_control_modules::common::TH1Reductor"};
    std::string moduleName{"QcCommon"};
  };

  /// Plotting options.
  struct PlotConfig {
    bool enablePerRingPlots{true};
    bool enableAllChannelsPlot{true};
    bool enableLegend{true};
    int  legendColumns{2};
  };

  std::vector<BeamConfig> beamConfigs;
  ChannelConfig           channelConfig;
  PlotConfig              plotConfig;

  /// Base colors per ring (ROOT color indices).
  std::vector<int> ringColors{44, 46, 38, 49, 29, 28};

  /// Marker styles within a ring (eight channels).
  std::vector<int> markerStyles{21, 22, 23, 24, 25, 26, 32, 33};

  /// Marker styles per ring (6 rings: 1..5.2) for the all-channels plot.
  std::vector<int> ringMarkerStyles{20, 21, 22, 23, 33, 34};
};

/// TrendingTaskConfig that generates a full FV0 setup from AutoTrendingConfig.
class AutoTrendingTaskConfig
  : public o2::quality_control::postprocessing::TrendingTaskConfig
{
 public:
  AutoTrendingTaskConfig() = default;
  AutoTrendingTaskConfig(std::string id, const boost::property_tree::ptree& config);
  ~AutoTrendingTaskConfig() = default;

  void generateFV0TrendingConfig(const AutoTrendingConfig& autoConfig);

 private:
  // Build data sources and plots
  void generateDataSources(const AutoTrendingConfig& autoConfig);
  void generatePerRingPlots(const AutoTrendingConfig& autoConfig,
                            const AutoTrendingConfig::BeamConfig& beam);
  void generateAllChannelsPlot(const AutoTrendingConfig& autoConfig,
                               const AutoTrendingConfig::BeamConfig& beam);

  // Small helpers
  int         getRingForChannel(int channel, int channelsPerRing) const;
  std::string getChannelName(int channel, int channelsPerRing) const;
  static std::string pad2(int v);
  static std::string yRangeFromGain(double g);

  // Graph/legend builders
  o2::quality_control::postprocessing::TrendingTaskConfig::Graph
  createChannelGraph(const std::string& beamName,
                     int channel,
                     const std::string& title,
                     const o2::quality_control::postprocessing::TrendingTaskConfig::GraphStyle& style) const;

  o2::quality_control::postprocessing::TrendingTaskConfig::LegendConfig
  createLegend(float x1, float y1, float x2, float y2, int columns, bool enabled = true) const;
};

} // namespace o2::quality_control_modules::fv0

#endif // QUALITYCONTROL_AUTOTRENDINGTASKCONFIG_H
