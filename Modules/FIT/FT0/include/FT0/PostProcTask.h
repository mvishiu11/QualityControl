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
/// \file   PostProcTask.h
/// \author Sebastian Bysiak sbysiak@cern.ch
///

#ifndef QC_MODULE_FT0_POSTPROCTASK_H
#define QC_MODULE_FT0_POSTPROCTASK_H

#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/DatabaseInterface.h"
#include "FITCommon/PostProcHelper.h"
#include "FITCommon/DetectorFIT.h"
#include "FT0/ChannelGeometry.h"

#include "CCDB/CcdbApi.h"
#include "CommonConstants/LHCConstants.h"
#include "FT0Base/Constants.h"
#include "DataFormatsFT0/ChannelData.h"
#include "DataFormatsFT0/Digit.h"

#include <TH2.h>
#include <TH2Poly.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <memory>
#include <string>

class TH1F;
class TCanvas;
class TLegend;
class TProfile;

namespace o2::quality_control_modules::ft0
{

/// \brief Basic Postprocessing Task for FT0, computes among others the trigger rates
/// \author Sebastian Bysiak sbysiak@cern.ch
class PostProcTask final : public quality_control::postprocessing::PostProcessingInterface
{
 public:
  PostProcTask() = default;
  ~PostProcTask() override;
  void configure(const boost::property_tree::ptree&) override;
  void initialize(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;
  void update(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;
  void finalize(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;
  constexpr static std::size_t sBCperOrbit = o2::constants::lhc::LHCMaxBunches;
  constexpr static std::size_t sNCHANNELS_PM = o2::ft0::Constants::sNCHANNELS_PM;
  using Detector_t = o2::quality_control_modules::fit::detectorFIT::DetectorFT0;

 private:
  o2::quality_control_modules::fit::PostProcHelper mPostProcHelper;
  bool mIsFirstIter{ true };
  typename Detector_t::TrgMap_t mMapPMbits = Detector_t::sMapPMbits;
  typename Detector_t::TrgMap_t mMapTechTrgBitsExtra = Detector_t::sMapTechTrgBitsExtra;
  typename Detector_t::TrgMap_t mMapTrgBits = Detector_t::sMapTrgBits;
  // MOs
  std::unique_ptr<TH2F> mHistChDataNOTbits;
  std::unique_ptr<TH1F> mHistTriggers;
  std::unique_ptr<TH1F> mHistTriggerRates;
  std::unique_ptr<TH1F> mHistTimeInWindow;
  std::unique_ptr<TH1F> mHistCFDEff;
  std::unique_ptr<TH1F> mHistChannelID_outOfBC;
  std::unique_ptr<TH1F> mHistTrg_outOfBC;
  std::unique_ptr<TH1F> mHistTrgValidation;
  std::unique_ptr<TH2F> mHistBcPattern;
  std::unique_ptr<TH2F> mHistBcTrgOutOfBunchColl;
  std::unique_ptr<TProfile> mAmpl;
  std::unique_ptr<TProfile> mTime;
  std::unique_ptr<TH2Poly> mHistStatsSideA;
  std::unique_ptr<TH2Poly> mHistStatsSideC;

  /// Sum of count in bins vs amplitude for detector channels 0 to 31 (A-side inner).
  /// In other words a projection of range 1-32 in the 2D histogram amplitude vs channel
  /// from the DigitQcTask.
  std::unique_ptr<TH1F> mHistAmpAInner;

  /// Sum of count in bins vs amplitude for detector channels 32 to 95 (A-side outer).
  /// In other words a projection of range 33-96 in the 2D histogram amplitude vs channel
  /// from the DigitQcTask.
  std::unique_ptr<TH1F> mHistAmpAOuter;

  /// Sum of count in bins vs amplitude for detector channels 96 to 207 (C-side).
  /// In other words a projection of range 97-208 in the 2D histogram amplitude vs channel
  /// from the DigitQcTask.
  std::unique_ptr<TH1F> mHistAmpC;

  /// Sum of normalized count in bins vs amplitude for all detector channels (0-207).
  /// In other words the sum of Y projections for each bin in the 2D histogram amplitude vs channel
  /// from the DigitQcTask, where each projection is scaled with 1/(counts in that projection)
  std::unique_ptr<TH1F> mHistAmpNormPerChannel;

  // NEW: MIP tracking histograms
  std::unique_ptr<TH1F> mHistMIPValues;     ///< MIP values per channel
  std::unique_ptr<TH1F> mHistMIPDeviations; ///< MIP deviations from expected
  std::unique_ptr<TH2F> mHistMIPTrends;     ///< MIP trends over time

  // NEW: Summary graphs for amplitude analysis
  std::unique_ptr<TGraphErrors> mGraphMPVvsChannel; ///< Mean vs Channel from Gaussian fits
  std::unique_ptr<TGraphErrors> mGraphMIPvsChannel; ///< MIP vs Channel summary

  ChannelGeometry mChannelGeometry; //!
  // Configurations
  int mLowTimeThreshold{ -192 };
  int mUpTimeThreshold{ 192 };
  std::string mAsynchChannelLogic{ "standard" };

  // NEW: Amplitude analysis configuration
  bool mEnableAmplitudeAnalysis{ false }; ///< Enable enhanced amplitude analysis
  bool mEnableMIPTracking{ false };       ///< Enable MIP tracking
  double mSliceFracLeft{ 0.25 };          ///< Left side fit window fraction
  double mSliceFracRight{ 0.25 };         ///< Right side fit window fraction
  double mMIPExpectedValue{ 100.0 };      ///< Expected MIP value in ADC channels
  double mMIPTolerancePercent{ 20.0 };    ///< Tolerance for MIP deviation (%)
  int mMinEntriesForFit{ 50 };            ///< Minimum entries required for Gaussian fit

  // NEW: Gaussian fit results storage
  std::array<double, sNCHANNELS_PM> mMean{};     ///< Fitted mean values
  std::array<double, sNCHANNELS_PM> mSigma{};    ///< Fitted sigma values
  std::array<double, sNCHANNELS_PM> mChanX{};    ///< Channel X coordinates for graphs
  std::array<double, sNCHANNELS_PM> mChanXErr{}; ///< Channel X error bars

  // NEW: MIP tracking arrays
  std::array<double, sNCHANNELS_PM> mMIPValues{};     ///< Current MIP values per channel
  std::array<double, sNCHANNELS_PM> mMIPDeviations{}; ///< MIP deviations from expected

  // NEW: Statistics tracking
  mutable int mSuccessfulFits{ 0 };        ///< Counter for successful fits
  mutable int mFailedFits{ 0 };            ///< Counter for failed fits
  mutable int mMIPChannelsInRange{ 0 };    ///< Counter for channels with MIP in expected range
  mutable int mMIPChannelsOutOfRange{ 0 }; ///< Counter for channels with MIP out of range

  void setTimestampToMOs();
  // TO REMOVE
  std::vector<unsigned int> mVecChannelIDs{};
  std::vector<std::string> mVecHistsToDecompose{};
  using HistDecomposed_t = TH1D;
  using MapHistsDecomposed_t = std::map<std::string, std::map<unsigned int, std::shared_ptr<HistDecomposed_t>>>;
  MapHistsDecomposed_t mMapHistsToDecompose{};
  void decomposeHists();
  void reset();

  // NEW: Amplitude analysis methods
  void performAmplitudeAnalysis(TH2F* hAmpPerChannel);
  void performMIPAnalysis();
  void updateMIPTrends();
  void logAmplitudeStatistics() const;
  std::pair<double, double> calculateFitWindow(double peak, double sliceFracLeft, double sliceFracRight) const;
};

} // namespace o2::quality_control_modules::ft0

#endif // QC_MODULE_FT0_POSTPROCTASK_H
