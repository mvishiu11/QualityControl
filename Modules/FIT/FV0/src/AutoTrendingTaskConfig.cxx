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
/// \file   AutoTrendingTaskConfig.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Lightweight knobs to auto-generate FV0 trending configuration.
///

#include "FV0/AutoTrendingTaskConfig.h"
#include <boost/property_tree/ptree.hpp>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace oqpp = o2::quality_control::postprocessing;
using TT = oqpp::TrendingTaskConfig;

namespace o2::quality_control_modules::fv0
{

std::string AutoTrendingTaskConfig::pad2(int v)
{
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(2) << v;
  return ss.str();
}

std::string AutoTrendingTaskConfig::yRangeFromGain(double g)
{
  const double lo = 0.5 * g;
  const double hi = 1.5 * g;
  std::ostringstream ss;
  ss.setf(std::ios::fixed);
  ss << std::setprecision(3) << lo << ":" << hi;
  return ss.str();
}

AutoTrendingTaskConfig::AutoTrendingTaskConfig(std::string id,
                                               const boost::property_tree::ptree& config)
  : TT(std::move(id), config)
{
  const auto key = "qc.postprocessing." + this->id + ".";

  const bool useAuto = config.get<bool>(key + "autoFV0Trending", false);
  if (!useAuto) {
    return;
  }

  AutoTrendingConfig autoCfg;

  // Behavior
  autoCfg.producePlotsOnUpdate = config.get<bool>(key + "producePlotsOnUpdate", true);
  autoCfg.resumeTrend          = config.get<bool>(key + "resumeTrend", true);
  autoCfg.trendIfAllInputs     = config.get<bool>(key + "trendIfAllInputs", false);
  autoCfg.trendingTimestamp    = config.get<std::string>(key + "trendingTimestamp", "validFrom");

  // Layout
  autoCfg.channelConfig.totalChannels   = config.get<int>(key + "totalChannels", 48);
  autoCfg.channelConfig.channelsPerRing = config.get<int>(key + "channelsPerRing", 8);
  autoCfg.channelConfig.dataSourcePath  = config.get<std::string>(
      key + "dataSourcePath", "FV0/MO/AmplitudePostProc/TrendsScalars");

  // Plot toggles
  autoCfg.plotConfig.enablePerRingPlots    = config.get<bool>(key + "enablePerRingPlots", true);
  autoCfg.plotConfig.enableAllChannelsPlot = config.get<bool>(key + "enableAllChannelsPlot", true);
  autoCfg.plotConfig.enableLegend          = config.get<bool>(key + "enableLegend", true);
  autoCfg.plotConfig.legendColumns         = config.get<int>(key + "legendColumns", 2);

  // Beam configs (expectedGain drives dynamic Y range)
  if (auto beams = config.get_child_optional(key + "beamConfigs")) {
    for (const auto& kv : beams.value()) {
      AutoTrendingConfig::BeamConfig b;
      b.name         = kv.second.get<std::string>("name");
      b.displayName  = kv.second.get<std::string>("displayName");
      b.expectedGain = kv.second.get<double>("expectedGain", 15.0);
      autoCfg.beamConfigs.push_back(std::move(b));
    }
  } else {
    AutoTrendingConfig::BeamConfig def;
    def.name         = "pp";
    def.displayName  = "pp Collisions";
    def.expectedGain = config.get<double>(key + "expectedGain", 15.0);
    autoCfg.beamConfigs.push_back(std::move(def));
  }

  // Generate full config
  generateFV0TrendingConfig(autoCfg);

  // Mirror to base
  producePlotsOnUpdate = autoCfg.producePlotsOnUpdate;
  resumeTrend          = autoCfg.resumeTrend;
  trendIfAllInputs     = autoCfg.trendIfAllInputs;
  trendingTimestamp    = autoCfg.trendingTimestamp;
}

void AutoTrendingTaskConfig::generateFV0TrendingConfig(const AutoTrendingConfig& autoConfig)
{
  plots.clear();
  dataSources.clear();

  generateDataSources(autoConfig);

  for (const auto& b : autoConfig.beamConfigs) {
    if (autoConfig.plotConfig.enablePerRingPlots) {
      generatePerRingPlots(autoConfig, b);
    }
    if (autoConfig.plotConfig.enableAllChannelsPlot) {
      generateAllChannelsPlot(autoConfig, b);
    }
  }
}

void AutoTrendingTaskConfig::generateDataSources(const AutoTrendingConfig& autoConfig)
{
  for (const auto& b : autoConfig.beamConfigs) {
    for (int ch = 0; ch < autoConfig.channelConfig.totalChannels; ++ch) {
      TT::DataSource src;
      src.type         = "repository";
      src.path         = autoConfig.channelConfig.dataSourcePath + "/" + b.name;
      src.name         = b.name + "_ch" + pad2(ch); // zero-padded
      src.reductorName = autoConfig.channelConfig.reductorName;
      src.moduleName   = autoConfig.channelConfig.moduleName;
      dataSources.push_back(std::move(src));
    }
  }
}

void AutoTrendingTaskConfig::generatePerRingPlots(const AutoTrendingConfig& ac,
                                                  const AutoTrendingConfig::BeamConfig& beam)
{
  const int rings = (ac.channelConfig.totalChannels + ac.channelConfig.channelsPerRing - 1) /
                    ac.channelConfig.channelsPerRing;

  for (int ring = 0; ring < rings; ++ring) {
    TT::Plot p;

    // Nested under beam in the CCDB viewer
    if (ring < 4) {
      p.name  = beam.name + "/ring" + std::to_string(ring + 1) + "_trend";
      p.title = "FV0 Ring " + std::to_string(ring + 1) + ": " + beam.displayName;
    } else {
      const int sub = ring - 3; // 1,2 for ring 5
      p.name  = beam.name + "/ring5" + std::to_string(sub) + "_trend";
      p.title = "FV0 Ring 5." + std::to_string(sub) + ": " + beam.displayName;
    }

    p.graphAxisLabel = "MIP ADC(ch):Time";
    p.graphYRange    = yRangeFromGain(beam.expectedGain);
    p.legend = createLegend(0.75f, 0.15f, 0.93f, 0.28f,
      std::max(1, ac.plotConfig.legendColumns), true);

    // Base ring color; channels get related shades via offsets.
    const int baseColor = (ring < static_cast<int>(ac.ringColors.size()))
                          ? ac.ringColors[ring] : 1;
    static const int kOffsets[8] = {0, +1, -1, +2, -2, +3, -3, +4};

    const int start = ring * ac.channelConfig.channelsPerRing;
    const int end   = std::min(start + ac.channelConfig.channelsPerRing,
                               ac.channelConfig.totalChannels);

    for (int ch = start; ch < end; ++ch) {
      TT::GraphStyle st;
      const int posInRing  = ch - start;
      const int colorShift = kOffsets[std::min(posInRing, 7)];
      const int lineColor  = std::max(1, baseColor + colorShift);

      st.lineColor   = lineColor;
      st.markerColor = lineColor;
      st.lineWidth   = 2;
      st.markerStyle = (posInRing < static_cast<int>(ac.markerStyles.size()))
                        ? ac.markerStyles[posInRing] : 20;

      const std::string chName = getChannelName(ch, ac.channelConfig.channelsPerRing);
      p.graphs.push_back(createChannelGraph(beam.name, ch, chName, st));
    }

    plots.push_back(std::move(p));
  }
}

void AutoTrendingTaskConfig::generateAllChannelsPlot(const AutoTrendingConfig& ac,
                                                     const AutoTrendingConfig::BeamConfig& beam)
{
  TT::Plot p;
  p.name           = beam.name + "/all_channels";
  p.title          = "FV0 All Channels: " + beam.displayName;
  p.graphAxisLabel = "MIP ADC(ch):Time";
  p.graphYRange    = yRangeFromGain(beam.expectedGain);
  p.legend = createLegend(0.15f, 0.15f, 0.90f, 0.30f, 8, true);

  const int rings = (ac.channelConfig.totalChannels + ac.channelConfig.channelsPerRing - 1) /
                    ac.channelConfig.channelsPerRing;

  for (int ring = 0; ring < rings; ++ring) {
    const int baseColor = (ring < static_cast<int>(ac.ringColors.size()))
                          ? ac.ringColors[ring] : 1;
    const int ringMarker = (ring < static_cast<int>(ac.ringMarkerStyles.size()))
                           ? ac.ringMarkerStyles[ring] : 20;

    const int start = ring * ac.channelConfig.channelsPerRing;
    const int end   = std::min(start + ac.channelConfig.channelsPerRing,
                               ac.channelConfig.totalChannels);

    for (int ch = start; ch < end; ++ch) {
      TT::GraphStyle st;
      st.lineColor   = baseColor;     // same color within ring
      st.markerColor = baseColor;
      st.markerStyle = ringMarker;    // same marker within ring
      st.lineWidth   = 1;

      // Legend/title shows A1–H1, A2–H2, … A52–H52
      const std::string title = getChannelName(ch, ac.channelConfig.channelsPerRing);
      p.graphs.push_back(createChannelGraph(beam.name, ch, title, st));
    }
  }

  plots.push_back(std::move(p));
}

int AutoTrendingTaskConfig::getRingForChannel(int channel, int channelsPerRing) const
{
  return channel / channelsPerRing;
}

std::string AutoTrendingTaskConfig::getChannelName(int channel, int channelsPerRing) const
{
  static const std::vector<std::string> sector{"A","B","C","D","E","F","G","H"};
  const int ring      = getRingForChannel(channel, channelsPerRing) + 1;
  const int posInRing = channel % channelsPerRing;

  if (ring <= 4) {
    return sector[posInRing] + std::to_string(ring);
  }
  // Ring 5.1 / 5.2
  const int sub = ring - 4;
  return sector[posInRing] + "5" + std::to_string(sub);
}

TT::Graph
AutoTrendingTaskConfig::createChannelGraph(const std::string& beamName,
                                           int channel,
                                           const std::string& title,
                                           const TT::GraphStyle& style) const
{
  TT::Graph g;
  const std::string cid = pad2(channel);

  g.name      = beamName + "_ch" + cid;
  g.title     = title;
  g.varexp    = beamName + "_ch" + cid + ".mean:time";
  g.selection = "";
  g.option    = "*LP";
  g.errors    = "";
  g.style     = style;

  return g;
}

TT::LegendConfig
AutoTrendingTaskConfig::createLegend(float x1, float y1, float x2, float y2, int columns, bool enabled) const
{
  TT::LegendConfig leg;
  leg.enabled  = enabled;
  leg.nColumns = columns;
  leg.x1 = x1; leg.y1 = y1;
  leg.x2 = x2; leg.y2 = y2;
  return leg;
}

} // namespace o2::quality_control_modules::fv0
