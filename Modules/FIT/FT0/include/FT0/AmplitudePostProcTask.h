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
/// \file   AmplitudePostProcTask.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Post-processing task for FT0 amplitude analysis
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
  // ---------- helpers ----------
  void reset();
  void setTimestampToMOs();

  // window and statistics
  std::pair<double, double> computeWindow(double peak) const;
  double computeWeightedMeanInWindow(const TH1D* h, double xmin, double xmax) const;

  // region fit (gaussian) with overlay on the region histogram
  bool fitRegionGaussianAndOverlay(TH1F* regionHist,
                                   const char* fitName,
                                   double& outMu,
                                   double& outSigma) const;

  // ---------- constants & config ----------
  static constexpr std::size_t sNCHANNELS_PM = o2::ft0::Constants::sNCHANNELS_PM;

  // histogram ADC range
  int mAmpMin{-100};
  int mAmpMax{4100};
  int mAmpBins{4200};

  // fixed fractional window (relative to peak position)
  double mLeftSliceFrac{0.15};
  double mRightSliceFrac{0.15};

  // weighted mean bin cut (ignore bins with < N entries)
  int mMinBinEntriesForWeight{1};

  // expected gain (ADC/MIP) for scaling / reference
  double mExpectedGain{14.0};

  // trending (persistent, per-run scalars)
  bool mTrendEnabled{true};
  int mTrendMaxPoints{10}; // (kept if you later add in-process graphs)
  std::string mTrendScalarsFolder{"TrendsScalars"}; // base folder inside our task path

  // ---------- managed objects ----------
  // per-channel histograms
  std::map<unsigned int, std::unique_ptr<TH1F>> mMapHistAmpPerChannel;

  // region histograms
  std::unique_ptr<TH1F> mHistAmpAll;     // 0-207
  std::unique_ptr<TH1F> mHistAmpAInner;  // 0-31
  std::unique_ptr<TH1F> mHistAmpAOuter;  // 32-95
  std::unique_ptr<TH1F> mHistAmpC;       // 96-207
  std::unique_ptr<TH1F> mHistAmpNormPerChannel;
  std::unique_ptr<TH1F> mHLastAInner, mHLastAOuter, mHLastC, mHLastAll;
  std::unique_ptr<TH1F> mTrendAInner, mTrendAOuter, mTrendC, mTrendAll;

  // per-channel summary (μW/expectedGain vs channel)
  std::unique_ptr<TGraphErrors> mGraphMeanNormVsChannel;

  // ---------- per-channel storage ----------
  std::array<double, sNCHANNELS_PM> mMeanW{};    // weighted mean
  std::array<double, sNCHANNELS_PM> mChanX{};
  std::array<double, sNCHANNELS_PM> mChanXErr{};

  // ---------- per-region last-fit results (Gaussian) ----------
  double mMuAInner{std::numeric_limits<double>::quiet_NaN()};
  double mMuAOuter{std::numeric_limits<double>::quiet_NaN()};
  double mMuC{std::numeric_limits<double>::quiet_NaN()};
  double mMuAll{std::numeric_limits<double>::quiet_NaN()};
  double mSigAInner{0.}, mSigAOuter{0.}, mSigC{0.}, mSigAll{0.};

  // path to DigitQcTask MOs (to fetch AmpPerChannel)
  std::string mPathDigitQcTask{"FT0/MO/Digits/"};

  // PostProc helper
  o2::quality_control_modules::fit::PostProcHelper mPostProcHelper;
};

} // namespace o2::quality_control_modules::ft0

#endif // QC_MODULE_FT0_AMPLITUDEPOSTPROCTASK_H
