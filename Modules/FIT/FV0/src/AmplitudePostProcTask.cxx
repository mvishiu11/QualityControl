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
/// \file   AmplitudePostProcTask.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Post-processing task for FV0 amplitude analysis per channel
///

// O2 QC / framework
#include "FV0/AmplitudePostProcTask.h"
#include "QualityControl/QcInfoLogger.h"
#include "FITCommon/HelperHist.h"
#include "FITCommon/HelperGraph.h"
#include "FITCommon/HelperCommon.h"
// ROOT
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
// STL
#include <limits>
#include <memory>
#include <sstream>

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::fit;

namespace o2::quality_control_modules::fv0
{

void AmplitudePostProcTask::configure(const boost::property_tree::ptree& cfg)
{
  const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
  const char* cfgCustom = Form("%s.custom", cfgPath);

  mPostProcHelper.configure(cfg, cfgPath, "FV0");

  auto cfgGet = [&cfgCustom](const std::string& entry) {
    return Form("%s.%s", cfgCustom, entry.c_str());
  };

  // Basic histogram configuration
  mAmpMin = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampMin"), -100);
  mAmpMax = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampMax"), 4100);
  mAmpBins = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampBins"), 4200);
  mSliceFrac = helper::getConfigFromPropertyTree<double>(cfg, cfgGet("sliceFrac"), 0.25);

  // Empirical fitting configuration
  mUseEmpiricalFitting = helper::getConfigFromPropertyTree<bool>(cfg, cfgGet("useEmpiricalFitting"), true);
  mUseFallbackFitting = helper::getConfigFromPropertyTree<bool>(cfg, cfgGet("useFallbackFitting"), true);

  ILOG(Info, Support) << "AmplitudePostProcTask configured with empirical fitting: "
                      << (mUseEmpiricalFitting ? "enabled" : "disabled")
                      << ", fallback: " << (mUseFallbackFitting ? "enabled" : "disabled")
                      << ", slice fraction: " << mSliceFrac << ENDM;
}

void AmplitudePostProcTask::reset()
{
  mMapHistAmpPerChannel.clear();
  mHistAmpAll.reset();
  mHistAmpNormPerChannel.reset();
  mGraphMPVDiv16.reset();

  mEmpiricalFitsUsed = 0;
  mFallbackFitsUsed = 0;
  mFailedFits = 0;
}

void AmplitudePostProcTask::setTimestampToMOs()
{
  for (int iObj = 0; iObj < getObjectsManager()->getNumberPublishedObjects(); iObj++) {
    auto mo = getObjectsManager()->getMonitorObject(iObj);
    mo->addOrUpdateMetadata(mPostProcHelper.mTimestampMetaField,
                            std::to_string(mPostProcHelper.mTimestampAnchor));
  }
}

void AmplitudePostProcTask::initializeEmpiricalParameters()
{
  ILOG(Info, Support) << "Initializing empirical slice fractions based on detector experience" << ENDM;

  // Initialize default slice fractions
  mDefaultFitParams = ChannelFitParams(0.28, 0.20, false, 1, "default_all_rings");

  // Create the channel mapping for detector position identification
  for (unsigned int ch = 0; ch < sNCHANNELS_PM && ch < 48; ++ch) {
    FV0Ring ring = static_cast<FV0Ring>(ch / 8);
    FV0Sector sector = static_cast<FV0Sector>(ch % 8);
    mChannelMapping[ch] = DetectorPosition(ring, sector);
  }

  // Ring R2 specific slice fractions
  // Channel 13 = Ring R2 (8-15), Sector F (+5) = 8+5 = 13
  mChannelFitParams[13] = ChannelFitParams(0.24, 0.28, false, 1, "F_R2");

  // Ring R4 specific slice fractions
  // Channel 29 = Ring R4 (24-31), Sector F (+5) = 24+5 = 29
  mChannelFitParams[29] = ChannelFitParams(0.24, 0.24, false, 1, "F_R4");
  // Channel 31 = Ring R4 (24-31), Sector H (+7) = 24+7 = 31
  mChannelFitParams[31] = ChannelFitParams(0.28, 0.25, false, 1, "H_R4");

  // Special case for all outer ring channels (32-47)
  for (unsigned int ch = 32; ch < 48; ++ch) {
    mChannelFitParams[ch] = ChannelFitParams(0.36, 0.36, false, 1, "Ring5x_default");
  }

  // Ring R51 specific overrides for channels with individual calibration
  mChannelFitParams[33] = ChannelFitParams(0.35, 0.45, false, 1, "B_R51"); // Ring R51, Sector B
  mChannelFitParams[34] = ChannelFitParams(0.26, 0.35, false, 1, "C_R51"); // Ring R51, Sector C
  mChannelFitParams[35] = ChannelFitParams(0.40, 0.45, false, 1, "D_R51"); // Ring R51, Sector D
  mChannelFitParams[36] = ChannelFitParams(0.36, 0.36, false, 1, "E_R51"); // Ring R51, Sector E
  mChannelFitParams[39] = ChannelFitParams(0.30, 0.40, false, 1, "H_R51"); // Ring R51, Sector H

  // Ring R52 specific overrides
  mChannelFitParams[40] = ChannelFitParams(0.30, 0.30, false, 1, "A_R52"); // Ring R52, Sector A
  mChannelFitParams[42] = ChannelFitParams(0.36, 0.36, false, 1, "C_R52"); // Ring R52, Sector C
  mChannelFitParams[43] = ChannelFitParams(0.36, 0.36, false, 1, "D_R52"); // Ring R52, Sector D
  mChannelFitParams[44] = ChannelFitParams(0.30, 0.30, false, 1, "E_R52"); // Ring R52, Sector E
  mChannelFitParams[45] = ChannelFitParams(0.36, 0.36, false, 1, "F_R52"); // Ring R52, Sector F
  mChannelFitParams[46] = ChannelFitParams(0.40, 0.40, false, 1, "G_R52"); // Ring R52, Sector G
  mChannelFitParams[47] = ChannelFitParams(0.34, 0.38, false, 1, "H_R52"); // Ring R52, Sector H

  ILOG(Info, Support) << "Initialized empirical slice fractions for "
                      << mChannelFitParams.size() << " channels with custom settings" << ENDM;
}

DetectorPosition AmplitudePostProcTask::getChannelPosition(unsigned int channel) const
{
  // Look up the detector position for this channel number
  auto it = mChannelMapping.find(channel);
  if (it != mChannelMapping.end()) {
    return it->second;
  }

  // Fallback calculation if channel not found in mapping
  FV0Ring ring = static_cast<FV0Ring>(std::min(channel / 8, 5U)); // Clamp to valid ring range
  FV0Sector sector = static_cast<FV0Sector>(channel % 8);

  ILOG(Debug, Support) << "Channel " << channel << " not found in mapping, using calculated position" << ENDM;
  return DetectorPosition(ring, sector);
}

ChannelFitParams AmplitudePostProcTask::getChannelFitParams(unsigned int channel) const
{
  // Check if we have specific empirical parameters for this channel
  auto it = mChannelFitParams.find(channel);
  if (it != mChannelFitParams.end()) {
    return it->second;
  }

  return mDefaultFitParams;
}

std::pair<double, double> AmplitudePostProcTask::calculateFitWindow(unsigned int channel,
                                                                    double peak, int peakBin,
                                                                    TH1D* histogram) const
{
  // Fallback to the original fractional window approach
  if (!mUseEmpiricalFitting) {
    double xmin = std::max<double>(peak - mSliceFrac * std::abs(peak), mAmpMin);
    double xmax = std::min<double>(peak + mSliceFrac * std::abs(peak), mAmpMax);
    mFallbackFitsUsed++;
    return std::make_pair(xmin, xmax);
  }

  // Get the empirical slice fractions for this specific channel
  ChannelFitParams params = getChannelFitParams(channel);
  DetectorPosition pos = getChannelPosition(channel);

  double xmin = std::max<double>(peak - params.leftSliceFrac * std::abs(peak),
                                 static_cast<double>(mAmpMin));
  double xmax = std::min<double>(peak + params.rightSliceFrac * std::abs(peak),
                                 static_cast<double>(mAmpMax));

  // Basic sanity check: ensure the window is mathematically valid
  if (xmax <= xmin) {
    if (mUseFallbackFitting) {
      ILOG(Warning, Support) << "Invalid empirical window for channel "
                             << channel << " (" << pos.toString()
                             << "), falling back to fractional window" << ENDM;

      double fallback_xmin = std::max<double>(peak - mSliceFrac * std::abs(peak), mAmpMin);
      double fallback_xmax = std::min<double>(peak + mSliceFrac * std::abs(peak), mAmpMax);
      mFallbackFitsUsed++;
      return std::make_pair(fallback_xmin, fallback_xmax);
    } else {
      ILOG(Error, Support) << "Invalid empirical window for channel " << channel
                           << ", but fallback is disabled" << ENDM;
      mFailedFits++;
    }
  }

  ILOG(Debug, Support) << "Channel " << channel << " (" << pos.toString()
                       << "): window [" << xmin << ", " << xmax
                       << "] using slice fractions [" << params.leftSliceFrac
                       << ", " << params.rightSliceFrac << "]" << ENDM;

  mEmpiricalFitsUsed++;
  return std::make_pair(xmin, xmax);
}

void AmplitudePostProcTask::logFittingStatistics() const
{
  int totalFits = mEmpiricalFitsUsed + mFallbackFitsUsed + mFailedFits;
  if (totalFits > 0) {
    ILOG(Info, Support) << "Fitting statistics: " << mEmpiricalFitsUsed << " empirical fits ("
                        << (100.0 * mEmpiricalFitsUsed / totalFits) << "%), "
                        << mFallbackFitsUsed << " fallback fits ("
                        << (100.0 * mFallbackFitsUsed / totalFits) << "%), "
                        << mFailedFits << " failed fits ("
                        << (100.0 * mFailedFits / totalFits) << "%)" << ENDM;
  }
}

void AmplitudePostProcTask::initialize(Trigger trig, framework::ServiceRegistryRef services)
{
  ILOG(Info, Support) << "Initialising AmplitudePostProcTask" << ENDM;
  mPostProcHelper.initialize(trig, services);
  reset();
  initializeEmpiricalParameters();

  // Create global histograms
  mHistAmpAll = helper::registerHist<TH1F>(
    getObjectsManager(),
    quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpAllChannels",
    "FV0 all channels;Channel amplitude (ADC ch);Counts",
    mAmpBins, mAmpMin, mAmpMax);

  mHistAmpNormPerChannel = helper::registerHist<TH1F>(
    getObjectsManager(),
    quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpNormPerChannel",
    "FV0 normalised (per channel);Channel amplitude (ADC ch);Normalised counts",
    mAmpBins, mAmpMin, mAmpMax);
  mHistAmpNormPerChannel->Sumw2(kFALSE);

  // Create per-channel histograms
  for (unsigned int ch = 0; ch < sNCHANNELS_PM; ++ch) {
    const std::string name = Form("AmplitudePerChannel/Amp_ch%02u", ch);
    const std::string title = Form("FV0 channel %u;Channel amplitude (ADC ch);Counts", ch);
    mMapHistAmpPerChannel[ch] = helper::registerHist<TH1F>(
      getObjectsManager(),
      quality_control::core::PublicationPolicy::ThroughStop,
      "", name, title,
      mAmpBins, mAmpMin, mAmpMax);
  }

  ILOG(Info, Support) << "AmplitudePostProcTask initialized with " << sNCHANNELS_PM
                      << " channels and empirical fitting for " << mChannelFitParams.size()
                      << " specially-tuned channels" << ENDM;
}

void AmplitudePostProcTask::update(Trigger trig, framework::ServiceRegistryRef serviceReg)
{
  mPostProcHelper.update(trig, serviceReg);

  constexpr double kScale = 15.0; // pp scaling factor

  auto h2 = mPostProcHelper.template getObject<TH2F>("AmpPerChannel");
  if (!h2) {
    ILOG(Error, Support) << "MO 'AmpPerChannel' not found in PostProcHelper" << ENDM;
    setTimestampToMOs();
    return;
  }

  ILOG(Debug, Support) << "Processing amplitude data with " << h2->GetEntries()
                       << " entries using " << (mUseEmpiricalFitting ? "empirical" : "standard")
                       << " fitting approach" << ENDM;

  mHistAmpAll->Reset();
  mHistAmpNormPerChannel->Reset();
  for (auto& [_, h] : mMapHistAmpPerChannel) {
    h->Reset();
  }

  mEmpiricalFitsUsed = 0;
  mFallbackFitsUsed = 0;
  mFailedFits = 0;

  for (int chBin = 1, nBins = h2->GetXaxis()->GetNbins(); chBin <= nBins; ++chBin) {
    unsigned int ch = chBin - 1;
    if (ch >= sNCHANNELS_PM)
      continue;

    std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("p_ch%02u", ch), chBin, chBin));
    proj->Sumw2(kFALSE);

    ChannelFitParams params = getChannelFitParams(ch);
    if (params.useRebin && params.rebinFactor > 1) {
      proj->Rebin(params.rebinFactor);
      ILOG(Debug, Support) << "Applied rebinning factor " << params.rebinFactor
                           << " to channel " << ch << ENDM;
    }

    mMapHistAmpPerChannel[ch]->Add(proj.get());
    mHistAmpAll->Add(proj.get());
    if (double intg = proj->Integral(); intg > 0.) {
      mHistAmpNormPerChannel->Add(proj.get(), 1.0 / intg);
    }

    // Perform Gaussian slice fit
    if (proj->GetEntries() > 50) {
      static TF1 fG("fG", "gaus", mAmpMin, mAmpMax);

      int bMax = proj->GetMaximumBin();
      double peak = proj->GetBinCenter(bMax);

      auto [xmin, xmax] = calculateFitWindow(ch, peak, bMax, proj.get());

      proj->Fit(&fG, "QNR", "", xmin, xmax);
      mMean[ch] = fG.GetParameter(1);
      mSigma[ch] = std::abs(fG.GetParameter(2));

      if (mMean[ch] < mAmpMin || mMean[ch] > mAmpMax || mSigma[ch] <= 0) {
        DetectorPosition pos = getChannelPosition(ch);
        ILOG(Warning, Support) << "Suspicious fit result for channel " << ch
                               << " (" << pos.toString() << "): mean=" << mMean[ch]
                               << ", sigma=" << mSigma[ch] << ENDM;
        mMean[ch] = std::numeric_limits<double>::quiet_NaN();
        mSigma[ch] = 0.;
        mFailedFits++;
      }
    } else {
      mMean[ch] = std::numeric_limits<double>::quiet_NaN();
      mSigma[ch] = 0.;
      ILOG(Debug, Support) << "Insufficient entries (" << proj->GetEntries()
                           << ") for channel " << ch << ", skipping fit" << ENDM;
    }

    mChanX[ch] = ch;
    mChanXErr[ch] = 0.0;
  }

  // Create or update the summary graph
  if (!mGraphMPVDiv16) {
    mGraphMPVDiv16 = helper::registerGraph<TGraphErrors>(
      getObjectsManager(),
      quality_control::core::PublicationPolicy::ThroughStop,
      "AP", "GaussianSummary/MeanVsChannel",
      "FV0: Scaled Gaussian mean;Channel ID;#mu / scaling factor (collision-dependent)",
      sNCHANNELS_PM);

    mGraphMPVDiv16->GetXaxis()->SetLimits(-1, 49);
    mGraphMPVDiv16->GetYaxis()->SetRangeUser(0., 2.);
    mGraphMPVDiv16->SetMarkerStyle(21);
    mGraphMPVDiv16->SetMarkerSize(1.1);
    mGraphMPVDiv16->SetMarkerColor(kRed);
    mGraphMPVDiv16->SetLineColor(kBlack);
    mGraphMPVDiv16->GetListOfFunctions()->Clear();
    auto* refLine = new TLine(-0.5, 1.0, 48.5, 1.0);
    refLine->SetLineColor(kBlue);
    refLine->SetLineStyle(2);
    refLine->SetLineWidth(2);
    mGraphMPVDiv16->GetListOfFunctions()->Add(refLine);

    ILOG(Info, Support) << "Created Gaussian summary graph with reference line" << ENDM;
  }

  for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
    mGraphMPVDiv16->SetPoint(i, mChanX[i], mMean[i] / kScale);
    mGraphMPVDiv16->SetPointError(i, mChanXErr[i], 0.);
  }

  logFittingStatistics();
  setTimestampToMOs();

  ILOG(Info, Support) << "AmplitudePostProcTask update completed successfully with empirical fitting" << ENDM;
}

void AmplitudePostProcTask::finalize(Trigger, framework::ServiceRegistryRef)
{
  int totalChannelsProcessed = 0;
  int totalEntriesProcessed = 0;

  for (auto& [ch, hist] : mMapHistAmpPerChannel) {
    if (hist && hist->GetEntries() > 0) {
      totalChannelsProcessed++;
      totalEntriesProcessed += hist->GetEntries();
      ILOG(Debug, Support) << "Channel " << ch << " processed "
                           << hist->GetEntries() << " entries" << ENDM;
    }
  }

  if (mUseEmpiricalFitting) {
    int customChannels = mChannelFitParams.size();
    int defaultChannels = totalChannelsProcessed - customChannels;

    ILOG(Info, Support) << "Empirical fitting summary: " << customChannels
                        << " channels with custom parameters, " << defaultChannels
                        << " channels using defaults" << ENDM;
  }

  logFittingStatistics();

  if (mGraphMPVDiv16) {
    ILOG(Info, Support) << "Gaussian summary graph contains "
                        << mGraphMPVDiv16->GetN() << " points" << ENDM;
  }

  ILOG(Info, Support) << "AmplitudePostProcTask finalized. Processed "
                      << totalChannelsProcessed << " channels with "
                      << totalEntriesProcessed << " total entries using "
                      << (mUseEmpiricalFitting ? "empirical" : "standard")
                      << " fitting approach." << ENDM;
}

} // namespace o2::quality_control_modules::fv0