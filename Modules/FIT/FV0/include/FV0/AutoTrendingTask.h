// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.
// GPLv3, see "COPYING".

#ifndef QUALITYCONTROL_AUTOTRENDINGTASK_H
#define QUALITYCONTROL_AUTOTRENDINGTASK_H

#include "QualityControl/TrendingTask.h"
#include "FV0/AutoTrendingTaskConfig.h"
#include <boost/property_tree/ptree_fwd.hpp>

namespace o2::quality_control_modules::fv0
{

/// A thin wrapper around TrendingTask that autogenerates config for FV0.
class AutoTrendingTask : public o2::quality_control::postprocessing::TrendingTask
{
 public:
  AutoTrendingTask() = default;
  ~AutoTrendingTask() override = default;

  void configure(const boost::property_tree::ptree& config) override;

 private:
  AutoTrendingTaskConfig mAutoConfig;
};

} // namespace o2::quality_control_modules::fv0

#endif // QUALITYCONTROL_AUTOTRENDINGTASK_H
