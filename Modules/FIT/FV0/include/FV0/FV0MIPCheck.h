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
/// \file   FV0MIPCheck.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  A simple check for FV0 MIP
///

#ifndef QC_MODULE_FV0_FV0MIPCHECK_H
#define QC_MODULE_FV0_FV0MIPCHECK_H

#include "QualityControl/CheckInterface.h"
#include "CommonConstants/LHCConstants.h"

namespace o2::quality_control_modules::fv0
{

class FV0MIPCheck : public o2::quality_control::checker::CheckInterface
{
 public:
  FV0MIPCheck()  = default;
  ~FV0MIPCheck() override = default;

  void configure() override;
  Quality check(std::map<std::string,
        std::shared_ptr<MonitorObject>>* moMap) override;
  void beautify(std::shared_ptr<MonitorObject> mo,
                Quality res = Quality::Null) override;
  std::string getAcceptedType() override { return "TGraphErrors"; }

 private:
  double mReferenceADC {15.};   ///< expected gain
  double mToleranceADC {1.};    ///< |y-ref| that triggers warning
  bool   mDrawWindow   {true};  ///< paint the band?

  double mMaxDeviation {-1.};   ///< filled in check(), shown in label
};

} // namespace o2::quality_control_modules::fv0

#endif
