// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General
// Public License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   AmplitudePostProcTask.h
/// \author Jakub Muszyński <jakub.milosz.muszynski@cern.ch>
///

#ifndef QC_MODULE_FT0_AMPLITUDEPOSTPROCTASK_H
#define QC_MODULE_FT0_AMPLITUDEPOSTPROCTASK_H

// O2 QC / framework
#include "QualityControl/PostProcessingInterface.h"
#include "FT0Base/Constants.h"
#include "FITCommon/PostProcHelper.h"

// ROOT
#include <TGraphErrors.h>
#include <TLine.h>
#include <TF1.h>

// STL
#include <array>
#include <map>
#include <memory>
#include <string>

class TH1F;
class TH2F;
class TH1D;

namespace o2::quality_control_modules::ft0
{

class AmplitudePostProcTask final : public quality_control::postprocessing::PostProcessingInterface
{
 public:
  AmplitudePostProcTask() = default;
  ~AmplitudePostProcTask() override = default;

  void configure(const boost::property_tree::ptree& config) override;
  void initialize(quality_control::postprocessing::Trigger trigger,
                  framework::ServiceRegistryRef services) override;
  void update(quality_control::postprocessing::Trigger trigger,
              framework::ServiceRegistryRef services) override;
  void finalize(quality_control::postprocessing::Trigger trigger,
                framework::ServiceRegistryRef services) override;

 private:
  void reset();
  void setTimestampToMOs();

  std::pair<double, double> computeWindow(double peak) const;
  double computeWeightedMeanInWindow(const TH1D* h, double xmin, double xmax) const;

  bool fitRegionGaussianAndOverlay(TH1F* regionHist,
                                   const char* fitName,
                                   double& outMu,
                                   double& outSigma) const;

  static constexpr std::size_t sNCHANNELS_PM = o2::ft0::Constants::sNCHANNELS_PM;

  // histogram ADC range
  int mAmpMin{ -100 };
  int mAmpMax{ 4100 };
  int mAmpBins{ 4200 };

  // fixed fractional window
  double mLeftSliceFrac{ 0.15 };
  double mRightSliceFrac{ 0.15 };

  // weighted mean bin cut
  int mMinBinEntriesForWeight{ 1 };

  // expected gain (ADC/MIP)
  double mExpectedGain{ 14.0 };

  // trending
  bool mTrendEnabled{ true };
  std::string mTrendScalarsFolder{ "TrendsScalars" };

  // ---------- managed objects ----------
  std::map<unsigned int, std::unique_ptr<TH1F>> mMapHistAmpPerChannel;

  std::unique_ptr<TH1F> mHistAmpAll;    // 0-207
  std::unique_ptr<TH1F> mHistAmpAInner; // 0-31
  std::unique_ptr<TH1F> mHistAmpAOuter; // 32-95
  std::unique_ptr<TH1F> mHistAmpC;      // 96-207
  std::unique_ptr<TH1F> mHistAmpNormPerChannel;
  std::unique_ptr<TH1F> mHLastAInner, mHLastAOuter, mHLastC, mHLastAll;
  std::unique_ptr<TH1F> mTrendAInner, mTrendAOuter, mTrendC, mTrendAll;

  std::unique_ptr<TGraphErrors> mGraphMeanNormVsChannel;

  std::array<double, sNCHANNELS_PM> mMeanW{};
  std::array<double, sNCHANNELS_PM> mChanX{};
  std::array<double, sNCHANNELS_PM> mChanXErr{};

  double mMuAInner{ std::numeric_limits<double>::quiet_NaN() };
  double mMuAOuter{ std::numeric_limits<double>::quiet_NaN() };
  double mMuC{ std::numeric_limits<double>::quiet_NaN() };
  double mMuAll{ std::numeric_limits<double>::quiet_NaN() };
  double mSigAInner{ 0. }, mSigAOuter{ 0. }, mSigC{ 0. }, mSigAll{ 0. };

  std::string mPathDigitQcTask{ "FT0/MO/Digits/" };

  o2::quality_control_modules::fit::PostProcHelper mPostProcHelper;
};

} // namespace o2::quality_control_modules::ft0

#endif // QC_MODULE_FT0_AMPLITUDEPOSTPROCTASK_H
