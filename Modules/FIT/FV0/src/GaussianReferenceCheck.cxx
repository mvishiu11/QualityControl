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
/// \file   GaussianReferenceCheck.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Checks if Gaussian mean values are consistent with reference line at 1.0
///

#include "FV0/GaussianReferenceCheck.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/Quality.h"
#include "QualityControl/QcInfoLogger.h"
#include "Common/Utils.h"

// ROOT includes
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TList.h>
#include <TMath.h>

#include <DataFormatsQualityControl/FlagType.h>
#include <DataFormatsQualityControl/FlagTypeFactory.h>

using namespace std;
using namespace o2::quality_control;

namespace o2::quality_control_modules::fv0
{

void GaussianReferenceCheck::configure()
{
  // Configure reference value
  if (auto param = mCustomParameters.find("referenceValue"); param != mCustomParameters.end()) {
    mReferenceValue = stof(param->second);
    ILOG(Debug, Support) << "configure() : using referenceValue = " << mReferenceValue << ENDM;
  } else {
    ILOG(Debug, Support) << "configure() : using default referenceValue = " << mReferenceValue << ENDM;
  }

  // Configure tolerance for warning level
  if (auto param = mCustomParameters.find("toleranceWarning"); param != mCustomParameters.end()) {
    mToleranceWarning = stof(param->second);
    ILOG(Debug, Support) << "configure() : using toleranceWarning = " << mToleranceWarning << ENDM;
  } else {
    ILOG(Debug, Support) << "configure() : using default toleranceWarning = " << mToleranceWarning << ENDM;
  }

  // Configure tolerance for error level
  if (auto param = mCustomParameters.find("toleranceError"); param != mCustomParameters.end()) {
    mToleranceError = stof(param->second);
    ILOG(Debug, Support) << "configure() : using toleranceError = " << mToleranceError << ENDM;
  } else {
    ILOG(Debug, Support) << "configure() : using default toleranceError = " << mToleranceError << ENDM;
  }

  // Configure minimum fraction of good channels
  if (auto param = mCustomParameters.find("minChannelFractionGood"); param != mCustomParameters.end()) {
    mMinChannelFractionGood = stof(param->second);
    ILOG(Debug, Support) << "configure() : using minChannelFractionGood = " << mMinChannelFractionGood << ENDM;
  } else {
    ILOG(Debug, Support) << "configure() : using default minChannelFractionGood = " << mMinChannelFractionGood << ENDM;
  }

  // Configure minimum fraction for medium quality
  if (auto param = mCustomParameters.find("minChannelFractionMedium"); param != mCustomParameters.end()) {
    mMinChannelFractionMedium = stof(param->second);
    ILOG(Debug, Support) << "configure() : using minChannelFractionMedium = " << mMinChannelFractionMedium << ENDM;
  } else {
    ILOG(Debug, Support) << "configure() : using default minChannelFractionMedium = " << mMinChannelFractionMedium << ENDM;
  }

  // Configure whether to require error bar overlap for good quality
  mRequireErrorBarOverlap = o2::quality_control_modules::common::getFromConfig(mCustomParameters, "requireErrorBarOverlap", true);
  mEnableMessage = o2::quality_control_modules::common::getFromConfig(mCustomParameters, "enableMessage", true);
}

Quality GaussianReferenceCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  Quality result = Quality::Null;
  TGraphErrors* graph = nullptr;

  // Find the Gaussian summary graph in the monitoring objects
  for (auto& [moName, mo] : *moMap) {
    if (mo->getName().find("MeanVsChannel") != std::string::npos ||
        mo->getName().find("GaussianSummary") != std::string::npos) {
      graph = dynamic_cast<TGraphErrors*>(mo->getObject());
      if (graph) {
        ILOG(Debug, Support) << "Found Gaussian summary graph: " << mo->getName() << ENDM;
        break;
      }
    }
  }

  // Check if we found the required graph
  if (!graph) {
    std::string reason = "Cannot find GaussianSummary/MeanVsChannel graph in monitoring objects";
    result.set(Quality::Null);
    result.addFlag(FlagTypeFactory::Unknown(), reason);
    ILOG(Warning) << reason << ENDM;
    return result;
  }

  // Initialize analysis variables
  mTotalChannelsAnalyzed = graph->GetN();
  mChannelsWithinReference = 0;
  mChannelsWithErrorBarOverlap = 0;
  mMeanDeviation = 0.0;
  mMaxDeviation = 0.0;

  // Analyze each point in the graph
  for (int i = 0; i < graph->GetN(); i++) {
    double x, y, ex, ey;
    graph->GetPoint(i, x, y);
    ex = graph->GetErrorX(i);
    ey = graph->GetErrorY(i);

    // Skip invalid points (NaN values from failed fits)
    if (TMath::IsNaN(y) || ey <= 0) {
      mTotalChannelsAnalyzed--; // Don't count invalid points
      continue;
    }

    // Calculate deviation from reference
    double deviation = TMath::Abs(y - mReferenceValue);
    mMeanDeviation += deviation;
    if (deviation > mMaxDeviation) {
      mMaxDeviation = deviation;
    }

    // Check if point is within tolerance (considering error bars)
    if (isPointWithinTolerance(y, ey, mToleranceWarning)) {
      mChannelsWithinReference++;
    }

    // Check if error bars overlap with reference line
    if (doesErrorBarOverlapReference(y, ey)) {
      mChannelsWithErrorBarOverlap++;
    }

    ILOG(Debug, Support) << "Channel " << i << ": value=" << y << " ± " << ey
                         << ", deviation=" << deviation
                         << ", overlap=" << doesErrorBarOverlapReference(y, ey) << ENDM;
  }

  // Calculate average deviation
  if (mTotalChannelsAnalyzed > 0) {
    mMeanDeviation /= mTotalChannelsAnalyzed;
  }

  // Determine quality based on analysis results
  result = Quality::Good;

  if (mTotalChannelsAnalyzed == 0) {
    result.set(Quality::Null);
    result.addFlag(FlagTypeFactory::Unknown(), "No valid channels found for analysis");
    return result;
  }

  // Calculate fractions for quality assessment
  float fractionWithinReference = static_cast<float>(mChannelsWithinReference) / mTotalChannelsAnalyzed;
  float fractionWithOverlap = static_cast<float>(mChannelsWithErrorBarOverlap) / mTotalChannelsAnalyzed;

  // Apply quality criteria
  if (fractionWithinReference < mMinChannelFractionMedium) {
    result.set(Quality::Bad);
    result.addFlag(FlagTypeFactory::Unknown(),
                   Form("Only %.1f%% of channels within tolerance (required: %.1f%% for medium quality)",
                        fractionWithinReference * 100, mMinChannelFractionMedium * 100));
  } else if (fractionWithinReference < mMinChannelFractionGood) {
    result.set(Quality::Medium);
    result.addFlag(FlagTypeFactory::Unknown(),
                   Form("Only %.1f%% of channels within tolerance (required: %.1f%% for good quality)",
                        fractionWithinReference * 100, mMinChannelFractionGood * 100));
  }

  // Additional check for error bar overlap if required
  if (mRequireErrorBarOverlap && result == Quality::Good) {
    if (fractionWithOverlap < mMinChannelFractionGood) {
      result.set(Quality::Medium);
      result.addFlag(FlagTypeFactory::Unknown(),
                     Form("Only %.1f%% of channels have error bars overlapping reference line",
                          fractionWithOverlap * 100));
    }
  }

  // Add metadata for trending and debugging
  result.addMetadata("totalChannelsAnalyzed", std::to_string(mTotalChannelsAnalyzed));
  result.addMetadata("channelsWithinReference", std::to_string(mChannelsWithinReference));
  result.addMetadata("channelsWithErrorBarOverlap", std::to_string(mChannelsWithErrorBarOverlap));
  result.addMetadata("meanDeviation", std::to_string(mMeanDeviation));
  result.addMetadata("maxDeviation", std::to_string(mMaxDeviation));

  ILOG(Info, Support) << "Gaussian reference check completed: "
                      << mChannelsWithinReference << "/" << mTotalChannelsAnalyzed
                      << " channels within tolerance, quality=" << result << ENDM;

  return result;
}

std::string GaussianReferenceCheck::getAcceptedType()
{
  return "TGraphErrors";
}

void GaussianReferenceCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
  auto* graph = dynamic_cast<TGraphErrors*>(mo->getObject());
  if (graph == nullptr) {
    ILOG(Warning, Devel) << "Could not cast " << mo->getName() << " to TGraphErrors*, skipping beautification" << ENDM;
    return;
  }

  if (!mEnableMessage) {
    return;
  }

  // Create informative message based on check results
  TPaveText* msg = new TPaveText(0.15, 0.02, 0.85, 0.15, "NDC");
  graph->GetListOfFunctions()->Add(msg);
  msg->SetName(Form("%s_msg", mo->GetName()));
  msg->Clear();
  msg->SetTextAlign(12);
  msg->SetTextSize(0.03);

  // Calculate summary statistics for display
  float fractionGood = mTotalChannelsAnalyzed > 0 ? static_cast<float>(mChannelsWithinReference) / mTotalChannelsAnalyzed : 0.0;
  float fractionOverlap = mTotalChannelsAnalyzed > 0 ? static_cast<float>(mChannelsWithErrorBarOverlap) / mTotalChannelsAnalyzed : 0.0;

  std::string summary = Form("Channels: %d/%d (%.1f%%) within tolerance, %d (%.1f%%) overlap reference line",
                             mChannelsWithinReference, mTotalChannelsAnalyzed, fractionGood * 100,
                             mChannelsWithErrorBarOverlap, fractionOverlap * 100);

  std::string deviation = Form("Mean deviation: %.3f, Max deviation: %.3f", mMeanDeviation, mMaxDeviation);

  // Set color and message based on quality
  if (checkResult == Quality::Good) {
    msg->SetFillColor(kGreen);
    msg->AddText((summary + " >> Quality::Good <<").c_str());
    msg->AddText(deviation.c_str());
  } else if (checkResult == Quality::Bad) {
    msg->SetFillColor(kRed);
    msg->AddText((summary + " >> Quality::Bad <<").c_str());
    msg->AddText(deviation.c_str());
    msg->AddText("Check detector calibration and hardware status!");
  } else if (checkResult == Quality::Medium) {
    msg->SetFillColor(kOrange);
    msg->AddText((summary + " >> Quality::Medium <<").c_str());
    msg->AddText(deviation.c_str());
    msg->AddText("Some channels showing deviations from reference");
  } else if (checkResult == Quality::Null) {
    msg->SetFillColor(kGray);
    msg->AddText("No valid data available for quality assessment >> Quality::Null <<");
  }
}

bool GaussianReferenceCheck::isPointWithinTolerance(double value, double error, double tolerance) const
{
  double effectiveTolerance = tolerance + error;
  return TMath::Abs(value - mReferenceValue) <= effectiveTolerance;
}

bool GaussianReferenceCheck::doesErrorBarOverlapReference(double value, double error) const
{
  return (value - error <= mReferenceValue) && (mReferenceValue <= value + error);
}

} // namespace o2::quality_control_modules::fv0