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
/// \brief  Post-processing task for FV0 amplitude analysis per channel
///

#ifndef QC_MODULE_FV0_AMPLITUDEPOSTPROCTASK_H
#define QC_MODULE_FV0_AMPLITUDEPOSTPROCTASK_H

// O2 QC / framework
#include "QualityControl/PostProcessingInterface.h"
#include "FV0Base/Constants.h"
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

namespace o2::quality_control_modules::fv0
{

// Empirical fitting parameters
struct ChannelFitParams {
  double leftSliceFrac;  ///< Left slice fraction (replaces mSliceFrac for left side)
  double rightSliceFrac; ///< Right slice fraction (replaces mSliceFrac for right side)
  bool useRebin;         ///< Whether to rebin histogram for this channel
  int rebinFactor;       ///< Rebin factor if rebinning is enabled
  std::string label;     ///< Human-readable label for debugging and logging

  /// Default constructor with fallback parameters (All rings params)
  ChannelFitParams()
    : leftSliceFrac(0.28), rightSliceFrac(0.20), useRebin(false), rebinFactor(1), label("default") {}

  /// Constructor with all parameters for custom channel configurations
  ChannelFitParams(double leftFrac, double rightFrac,
                   bool rebin = false, int rebinFac = 1, const std::string& lbl = "custom")
    : leftSliceFrac(leftFrac), rightSliceFrac(rightFrac), useRebin(rebin), rebinFactor(rebinFac), label(lbl) {}
};

/// Detector geometry mapping for FV0
/// These enums help translate between linear channel numbers and physical detector positions
enum class FV0Ring { R1 = 0,
                     R2 = 1,
                     R3 = 2,
                     R4 = 3,
                     R51 = 4,
                     R52 = 5 };
enum class FV0Sector { A = 0,
                       B = 1,
                       C = 2,
                       D = 3,
                       E = 4,
                       F = 5,
                       G = 6,
                       H = 7 };

/// Helper structure to represent detector position
struct DetectorPosition {
  FV0Ring ring;
  FV0Sector sector;

  /// Default constructor required for ROOT serialization
  /// Initializes to Ring R1, Sector A as a safe default
  DetectorPosition() : ring(FV0Ring::R1), sector(FV0Sector::A) {}

  /// Constructor with explicit ring and sector specification
  DetectorPosition(FV0Ring r, FV0Sector s) : ring(r), sector(s) {}

  /// Generate human-readable string for debugging and logging
  /// Example output: "D_R51" for Ring 51, Sector D
  std::string toString() const
  {
    const char* ringNames[] = { "R1", "R2", "R3", "R4", "R51", "R52" };
    const char* sectorNames[] = { "A", "B", "C", "D", "E", "F", "G", "H" };
    return std::string(sectorNames[static_cast<int>(sector)]) + "_" +
           std::string(ringNames[static_cast<int>(ring)]);
  }
};

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
  o2::quality_control_modules::fit::PostProcHelper mPostProcHelper;

  void reset();
  void setTimestampToMOs();

  // Empirical fitting
  void initializeEmpiricalParameters();
  DetectorPosition getChannelPosition(unsigned int channel) const;
  ChannelFitParams getChannelFitParams(unsigned int channel) const;
  std::pair<double, double> calculateFitWindow(unsigned int channel, double peak,
                                               int peakBin, TH1D* histogram) const;
  void logFittingStatistics() const;

  // Configuration & constants
  static constexpr std::size_t sNCHANNELS_PM = o2::fv0::Constants::nFv0ChannelsPlusRef;

  // Configuration parameters
  int mAmpMin{ -100 };               ///< Histogram ADC min
  int mAmpMax{ 4100 };               ///< Histogram ADC max
  int mAmpBins{ 4200 };              ///< Histogram bins
  double mSliceFrac{ 0.25 };         ///< Fallback fit window half-width (fraction of peak)
  bool mUseEmpiricalFitting{ true }; ///< Enable empirical fitting parameters
  bool mUseFallbackFitting{ true };  ///< Enable fallback to fractional window

  // QC managed objects
  std::map<unsigned int, std::unique_ptr<TH1F>> mMapHistAmpPerChannel;
  std::unique_ptr<TH1F> mHistAmpAll;
  std::unique_ptr<TH1F> mHistAmpNormPerChannel;

  // Gaussian fit results
  std::array<double, sNCHANNELS_PM> mMean{};
  std::array<double, sNCHANNELS_PM> mSigma{};
  std::array<double, sNCHANNELS_PM> mChanX{};
  std::array<double, sNCHANNELS_PM> mChanXErr{};

  // Gaussian error graph
  std::unique_ptr<TGraphErrors> mGraphMPVDiv16;

  // Empirical fitting data structures
  std::map<unsigned int, ChannelFitParams> mChannelFitParams; ///< Channel-specific fitting parameters
  ChannelFitParams mDefaultFitParams;                         ///< Default parameters for unmapped channels
  std::map<unsigned int, DetectorPosition> mChannelMapping;   ///< Maps channel numbers to detector positions

  // Statistics tracking
  mutable int mEmpiricalFitsUsed{ 0 }; ///< Counter for empirical fits applied
  mutable int mFallbackFitsUsed{ 0 };  ///< Counter for fallback fits applied
  mutable int mFailedFits{ 0 };        ///< Counter for failed fits
};

} // namespace o2::quality_control_modules::fv0

#endif // QC_MODULE_FV0_AMPLITUDEPOSTPROCTASK_H