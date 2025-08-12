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
/// \file   AutoTrendingTask.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Thin wrapper around TrendingTask to auto-generate FV0 trending configuration.
///

#ifndef QUALITYCONTROL_AUTOTRENDINGTASK_H
#define QUALITYCONTROL_AUTOTRENDINGTASK_H

#include "QualityControl/TrendingTask.h"
#include "FV0/AutoTrendingTaskConfig.h"
#include <boost/property_tree/ptree_fwd.hpp>

namespace o2::quality_control_modules::fv0
{

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
