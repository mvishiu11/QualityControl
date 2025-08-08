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
/// \file   AutoTrendingTask.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Thin wrapper around TrendingTask to auto-generate FV0 trending configuration.
///

#include "FV0/AutoTrendingTask.h"
#include "FV0/AutoTrendingTaskConfig.h"
#include "QualityControl/QcInfoLogger.h"

#include <boost/property_tree/ptree.hpp>

using o2::quality_control::postprocessing::TrendingTask;

namespace o2::quality_control_modules::fv0
{

// Serialize mAutoConfig back into the ptree schema TrendingTaskConfig expects.
static void serializeConfigToPtree(const AutoTrendingTaskConfig& cfg,
                                   boost::property_tree::ptree& root,
                                   const std::string& id)
{
  const std::string base = "qc.postprocessing." + id + ".";

  root.put(base + "producePlotsOnUpdate", cfg.producePlotsOnUpdate);
  root.put(base + "resumeTrend", cfg.resumeTrend);
  root.put(base + "trendIfAllInputs", cfg.trendIfAllInputs);
  root.put(base + "trendingTimestamp", cfg.trendingTimestamp);

  // dataSources[]
  boost::property_tree::ptree dsArray;
  for (const auto& ds : cfg.dataSources) {
    boost::property_tree::ptree n;
    n.put("type", ds.type);
    n.put("path", ds.path);
    n.put("name", ds.name);
    n.put("reductorName", ds.reductorName);
    n.put("moduleName", ds.moduleName);
    // Note: ds.reductorParameters is a CustomParameters object; omit here -> defaults.
    dsArray.push_back(std::make_pair("", n));
  }
  root.put_child(base + "dataSources", dsArray);

  // plots[]
  boost::property_tree::ptree plotsArray;
  for (const auto& p : cfg.plots) {
    boost::property_tree::ptree pn;
    pn.put("name", p.name);
    pn.put("title", p.title);
    pn.put("graphAxisLabel", p.graphAxisLabel);
    pn.put("graphYRange", p.graphYRange);
    pn.put("colorPalette", p.colorPalette);

    // legend
    boost::property_tree::ptree leg;
    leg.put("enabled", p.legend.enabled);
    leg.put("nColumns", p.legend.nColumns);
    leg.put("x1", p.legend.x1);
    leg.put("y1", p.legend.y1);
    leg.put("x2", p.legend.x2);
    leg.put("y2", p.legend.y2);
    pn.add_child("legend", leg);

    // graphs[]
    boost::property_tree::ptree graphsArray;
    for (const auto& g : p.graphs) {
      boost::property_tree::ptree gn;
      gn.put("name", g.name);
      gn.put("title", g.title);
      gn.put("varexp", g.varexp);
      gn.put("selection", g.selection);
      gn.put("option", g.option);
      gn.put("errors", g.errors);

      // style
      boost::property_tree::ptree st;
      st.put("lineColor", g.style.lineColor);
      st.put("lineStyle", g.style.lineStyle);
      st.put("lineWidth", g.style.lineWidth);
      st.put("markerColor", g.style.markerColor);
      st.put("markerStyle", g.style.markerStyle);
      st.put("markerSize", g.style.markerSize);
      st.put("fillColor", g.style.fillColor);
      st.put("fillStyle", g.style.fillStyle);
      gn.add_child("style", st);

      graphsArray.push_back(std::make_pair("", gn));
    }
    pn.add_child("graphs", graphsArray);

    plotsArray.push_back(std::make_pair("", pn));
  }
  root.put_child(base + "plots", plotsArray);
}

void AutoTrendingTask::configure(const boost::property_tree::ptree& config)
{
  // Build auto config from the incoming ptree
  mAutoConfig = AutoTrendingTaskConfig(getID(), config);

  // Re-emit a full config tree for the base class to parse
  auto fullCfg = config; // keep user keys intact
  serializeConfigToPtree(mAutoConfig, fullCfg, getID());

  ILOG(Info, Support) << "AutoTrendingTask: generated "
                      << mAutoConfig.dataSources.size() << " dataSources and "
                      << mAutoConfig.plots.size() << " plots; delegating to TrendingTask."
                      << ENDM;

  // Let the base do its normal setup (fills mConfig, builds reductors, etc.)
  TrendingTask::configure(fullCfg);
}

} // namespace o2::quality_control_modules::fv0
