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
/// \file   PostProcTask.cxx
/// \author Artur Furs afurs@cern.ch
///

#include "FT0/PostProcTask.h"
#include "QualityControl/QcInfoLogger.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsFT0/ChannelData.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TF1.h>
#include <TLine.h>
#include <TMath.h>

#include "FITCommon/HelperHist.h"
#include "FITCommon/HelperGraph.h"
#include "FITCommon/HelperCommon.h"

#include <iostream>
#include <string>
#include <limits>

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::fit;

namespace o2::quality_control_modules::ft0
{

PostProcTask::~PostProcTask()
{
}

void PostProcTask::configure(const boost::property_tree::ptree& config)
{
  const char* configPath = Form("qc.postprocessing.%s", getID().c_str());
  const char* configCustom = Form("%s.custom", configPath);
  auto cfgPath = [&configCustom](const std::string& entry) {
    return Form("%s.%s", configCustom, entry.c_str());
  };
  mPostProcHelper.configure(config, configPath, "FT0");

  // Detector related configs
  mLowTimeThreshold = helper::getConfigFromPropertyTree<int>(config, cfgPath("lowTimeThreshold"), -192);
  mUpTimeThreshold = helper::getConfigFromPropertyTree<int>(config, cfgPath("upTimeThreshold"), 192);
  mAsynchChannelLogic = helper::getConfigFromPropertyTree<std::string>(config, cfgPath("asynchChannelLogic"), "standard");
  mIsFirstIter = true; // to be sure

  // NEW: Amplitude analysis configuration
  mEnableAmplitudeAnalysis = helper::getConfigFromPropertyTree<bool>(config, cfgPath("enableAmplitudeAnalysis"), false);
  mEnableMIPTracking = helper::getConfigFromPropertyTree<bool>(config, cfgPath("enableMIPTracking"), false);
  mSliceFracLeft = helper::getConfigFromPropertyTree<double>(config, cfgPath("sliceFracLeft"), 0.25);
  mSliceFracRight = helper::getConfigFromPropertyTree<double>(config, cfgPath("sliceFracRight"), 0.25);
  mMIPExpectedValue = helper::getConfigFromPropertyTree<double>(config, cfgPath("mipExpectedValue"), 14.0);
  mMIPTolerancePercent = helper::getConfigFromPropertyTree<double>(config, cfgPath("mipTolerancePercent"), 20.0);
  mMinEntriesForFit = helper::getConfigFromPropertyTree<int>(config, cfgPath("minEntriesForFit"), 50);

  ILOG(Info, Support) << "PostProcTask configured with amplitude analysis: " 
                      << (mEnableAmplitudeAnalysis ? "enabled" : "disabled") 
                      << ", MIP tracking: " << (mEnableMIPTracking ? "enabled" : "disabled")
                      << ", slice fractions: [" << mSliceFracLeft << ", " << mSliceFracRight << "]" << ENDM;

  // TO REMOVE
  // VERY BAD SOLUTION, YOU SHOULDN'T USE IT
  const std::string del = ",";
  const std::string strChannelIDs = helper::getConfigFromPropertyTree<std::string>(config, cfgPath("channelIDs"), "");
  const std::string strHistsToDecompose = helper::getConfigFromPropertyTree<std::string>(config, cfgPath("histsToDecompose"), "");
  if (strChannelIDs.size() > 0 && strHistsToDecompose.size() > 0) {
    mVecChannelIDs = helper::parseParameters<unsigned int>(strChannelIDs, del);
    mVecHistsToDecompose = helper::parseParameters<std::string>(strHistsToDecompose, del);
  }
}

void PostProcTask::reset()
{
  mHistChDataNOTbits.reset();
  mHistTriggers.reset();
  mHistTriggerRates.reset();
  mHistTimeInWindow.reset();
  mHistCFDEff.reset();
  mHistChannelID_outOfBC.reset();
  mHistTrg_outOfBC.reset();
  mHistTrgValidation.reset();
  mHistBcPattern.reset();
  mHistBcTrgOutOfBunchColl.reset();
  mAmpl.reset();
  mTime.reset();
  mHistAmpAInner.reset();
  mHistAmpAOuter.reset();
  mHistAmpC.reset();
  mHistAmpNormPerChannel.reset();
  mMapHistsToDecompose.clear();

  // NEW: Reset amplitude analysis objects
  if (mEnableAmplitudeAnalysis) {
    mGraphMPVvsChannel.reset();
  }
  if (mEnableMIPTracking) {
    mHistMIPValues.reset();
    mHistMIPDeviations.reset();
    mHistMIPTrends.reset();
    mGraphMIPvsChannel.reset();
  }
  
  // Reset statistics counters
  mSuccessfulFits = 0;
  mFailedFits = 0;
  mMIPChannelsInRange = 0;
  mMIPChannelsOutOfRange = 0;
}

void PostProcTask::initialize(Trigger trg, framework::ServiceRegistryRef services)
{
  reset(); // for SSS procedure
  mPostProcHelper.initialize(trg, services);
  mIsFirstIter = true; // to be sure

  mHistChDataNOTbits = helper::registerHist<TH2F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "COLZ", "ChannelDataNegBits", "ChannelData NOT PM bits per ChannelID;Channel;Negative bit", sNCHANNELS_PM, 0, sNCHANNELS_PM, mMapPMbits);
  mHistTriggers = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "Triggers", "Triggers from TCM", mMapTechTrgBitsExtra);
  mHistTriggerRates = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "HIST", "TriggerRates", "Trigger rates; Triggers; Rate [kHz]", mMapTechTrgBitsExtra);
  mHistBcTrgOutOfBunchColl = helper::registerHist<TH2F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "COLZ", "OutOfBunchColl_BCvsTrg", "BC vs Triggers for out-of-bunch collisions;BC;Triggers", sBCperOrbit, 0, sBCperOrbit, mMapTechTrgBitsExtra);
  mHistBcPattern = helper::registerHist<TH2F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "COLZ", "bcPattern", "BC pattern", sBCperOrbit, 0, sBCperOrbit, mMapTechTrgBitsExtra);
  mHistTimeInWindow = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "TimeInWindowFraction", Form("Fraction of events with CFD in time gate(%i,%i) vs ChannelID;ChannelID;Event fraction with CFD in time gate", mLowTimeThreshold, mUpTimeThreshold), sNCHANNELS_PM, 0, sNCHANNELS_PM);
  mHistCFDEff = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "CFD_efficiency", "Fraction of events with CFD in ADC gate vs ChannelID;ChannelID;Event fraction with CFD in ADC gate;", sNCHANNELS_PM, 0, sNCHANNELS_PM);
  mHistChannelID_outOfBC = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "ChannelID_outOfBC", "ChannelID, out of bunch", sNCHANNELS_PM, 0, sNCHANNELS_PM);
  mHistTrg_outOfBC = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "Triggers_outOfBC", "Trigger fraction, out of bunch", mMapTechTrgBitsExtra);
  mHistTrgValidation = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "TrgValidation", "SW + HW only to validated triggers fraction", mMapTrgBits);
  mAmpl = helper::registerHist<TProfile>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "MeanAmplPerChannel", "mean ampl per channel;Channel;Ampl #mu #pm #sigma", o2::ft0::Constants::sNCHANNELS_PM, 0, o2::ft0::Constants::sNCHANNELS_PM);
  mTime = helper::registerHist<TProfile>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "MeanTimePerChannel", "mean time per channel;Channel;Time #mu #pm #sigma", o2::ft0::Constants::sNCHANNELS_PM, 0, o2::ft0::Constants::sNCHANNELS_PM);
  mHistAmpAInner = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "AmpAInner", "FT0 A-side inner channel amplitudes (ch 0-31);Channel amplitude (ADC ch);Counts", 4200, -100, 4100);
  mHistAmpAOuter = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "AmpAOuter", "FT0 A-side outer channel amplitudes (ch 32-95);Channel amplitude (ADC ch);Counts", 4200, -100, 4100);
  mHistAmpC = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "AmpC", "FT0 C-side channel ampltitudes (ch 96-207);Channel amplitude (ADC ch);Counts", 4200, -100, 4100);
  mHistAmpNormPerChannel = helper::registerHist<TH1F>(getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop, "", "AmpNormPerChannel", "FT0 channel amplitudes normalized per channel (ch 0-207);Channel amplitude (ADC ch);#sum_{ch}AmpHist_{ch} #times (1/counts_{ch})", 4200, -100, 4100);
  mHistAmpNormPerChannel->Sumw2(kFALSE);

  // NEW: Create MIP tracking histograms if enabled
  if (mEnableMIPTracking) {
    mHistMIPValues = helper::registerHist<TH1F>(
        getObjectsManager(),
        quality_control::core::PublicationPolicy::ThroughStop,
        "", "MIPValues",
        "MIP values per channel;Channel;MIP value (ADC ch)",
        sNCHANNELS_PM, 0, sNCHANNELS_PM);
        
    mHistMIPDeviations = helper::registerHist<TH1F>(
        getObjectsManager(),
        quality_control::core::PublicationPolicy::ThroughStop,
        "", "MIPDeviations",
        "MIP deviations from expected;Channel;Deviation (%)",
        sNCHANNELS_PM, 0, sNCHANNELS_PM);
        
    mHistMIPTrends = helper::registerHist<TH2F>(
        getObjectsManager(),
        quality_control::core::PublicationPolicy::ThroughStop,
        "COLZ", "MIPTrends",
        "MIP trends over time;Time bin;Channel;MIP deviation (%)",
        100, 0, 100, sNCHANNELS_PM, 0, sNCHANNELS_PM);
        
    ILOG(Info, Support) << "Created MIP tracking histograms" << ENDM;
  }

  mChannelGeometry.init(-200., 200., -200., 200., 10.); // values - borders for hist and margin

  mHistStatsSideA = mChannelGeometry.makeHistSideA("GeoChannelStatA", "Channel occupancy, side-A");
  mHistStatsSideC = mChannelGeometry.makeHistSideC("GeoChannelStatC", "Channel occupancy, side-C");
  getObjectsManager()->startPublishing(mHistStatsSideA.get());
  getObjectsManager()->setDefaultDrawOptions(mHistStatsSideA.get(), "TEXT COLZ L");
  getObjectsManager()->startPublishing(mHistStatsSideC.get());
  getObjectsManager()->setDefaultDrawOptions(mHistStatsSideC.get(), "TEXT COLZ L");

  ILOG(Info, Support) << "PostProcTask initialized with " << sNCHANNELS_PM 
                      << " channels, amplitude analysis " << (mEnableAmplitudeAnalysis ? "enabled" : "disabled")
                      << ", MIP tracking " << (mEnableMIPTracking ? "enabled" : "disabled") << ENDM;
}

std::pair<double, double> PostProcTask::calculateFitWindow(double peak, double sliceFracLeft, double sliceFracRight) const
{
    // Calculate asymmetric fit window using configurable fractions
    double xmin = std::max<double>(peak - sliceFracLeft * std::abs(peak), -100.0);
    double xmax = std::min<double>(peak + sliceFracRight * std::abs(peak), 4100.0);
    
    return std::make_pair(xmin, xmax);
}

void PostProcTask::performAmplitudeAnalysis(TH2F* hAmpPerChannel)
{
    if (!hAmpPerChannel || !mEnableAmplitudeAnalysis) return;
    
    ILOG(Debug, Support) << "Performing amplitude analysis with " << hAmpPerChannel->GetEntries() 
                         << " entries, slice fractions [" << mSliceFracLeft << ", " << mSliceFracRight << "]" << ENDM;
    
    // Reset statistics counters
    mSuccessfulFits = 0;
    mFailedFits = 0;
    
    // Process each channel with Gaussian fitting
    for (int chBin = 1, nBins = hAmpPerChannel->GetXaxis()->GetNbins(); chBin <= nBins; ++chBin) {
        unsigned int ch = chBin - 1;
        if (ch >= sNCHANNELS_PM) continue;
        
        // Extract amplitude histogram for this channel
        std::unique_ptr<TH1D> proj(hAmpPerChannel->ProjectionY(Form("p_ch%03u", ch), chBin, chBin));
        proj->Sumw2(kFALSE);
        
        // Initialize arrays
        mChanX[ch] = ch;
        mChanXErr[ch] = 0.0;
        
        // Perform Gaussian slice fit if we have enough statistics
        if (proj->GetEntries() > mMinEntriesForFit) {
            static TF1 fG("fG", "gaus", -100, 4100);
            
            int bMax = proj->GetMaximumBin();
            double peak = proj->GetBinCenter(bMax);
            
            // Use configurable fit window calculation
            auto [xmin, xmax] = calculateFitWindow(peak, mSliceFracLeft, mSliceFracRight);
            
            // Perform the Gaussian fit
            proj->Fit(&fG, "QNR", "", xmin, xmax);
            mMean[ch] = fG.GetParameter(1);
            mSigma[ch] = std::abs(fG.GetParameter(2));
            
            // Validate fit results
            if (mMean[ch] < -100 || mMean[ch] > 4100 || mSigma[ch] <= 0) {
                ILOG(Debug, Support) << "Invalid fit result for channel " << ch 
                                     << ": mean=" << mMean[ch] << ", sigma=" << mSigma[ch] << ENDM;
                mMean[ch] = std::numeric_limits<double>::quiet_NaN();
                mSigma[ch] = 0.;
                mFailedFits++;
            } else {
                mSuccessfulFits++;
            }
        } else {
            mMean[ch] = std::numeric_limits<double>::quiet_NaN();
            mSigma[ch] = 0.;
            ILOG(Debug, Support) << "Insufficient entries (" << proj->GetEntries() 
                                 << ") for channel " << ch << ", skipping fit" << ENDM;
        }
    }
    
    // Create or update summary graph
    if (!mGraphMPVvsChannel) {
        mGraphMPVvsChannel = helper::registerGraph<TGraphErrors>(
            getObjectsManager(),
            quality_control::core::PublicationPolicy::ThroughStop,
            "AP", "AmplitudeAnalysis/MeanVsChannel",
            "FT0: Gaussian mean per channel;Channel ID;Mean amplitude (ADC ch)",
            sNCHANNELS_PM);
            
        mGraphMPVvsChannel->GetXaxis()->SetLimits(-1, sNCHANNELS_PM + 1);
        mGraphMPVvsChannel->SetMarkerStyle(21);
        mGraphMPVvsChannel->SetMarkerSize(1.1);
        mGraphMPVvsChannel->SetMarkerColor(kRed);
        mGraphMPVvsChannel->SetLineColor(kBlack);
        
        ILOG(Info, Support) << "Created amplitude analysis summary graph" << ENDM;
    }
    
    // Update the graph data
    for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
        mGraphMPVvsChannel->SetPoint(i, mChanX[i], mMean[i]);
        mGraphMPVvsChannel->SetPointError(i, mChanXErr[i], 0.);
    }
}

void PostProcTask::performMIPAnalysis()
{
    if (!mEnableMIPTracking) return;
    
    mMIPChannelsInRange = 0;
    mMIPChannelsOutOfRange = 0;
    
    for (unsigned int ch = 0; ch < sNCHANNELS_PM; ++ch) {
        // Use the fitted mean as our MIP estimate
        mMIPValues[ch] = mMean[ch];
        
        // Calculate deviation from expected MIP value
        if (!std::isnan(mMean[ch]) && mMean[ch] > 0) {
            double deviation = ((mMean[ch] - mMIPExpectedValue) / mMIPExpectedValue) * 100.0;
            mMIPDeviations[ch] = deviation;
            
            // Check if MIP is within expected range
            if (std::abs(deviation) <= mMIPTolerancePercent) {
                mMIPChannelsInRange++;
            } else {
                mMIPChannelsOutOfRange++;
                ILOG(Debug, Support) << "Channel " << ch << " MIP deviation: " << deviation 
                                     << "% (mean=" << mMean[ch] << ", expected=" << mMIPExpectedValue << ")" << ENDM;
            }
        } else {
            mMIPDeviations[ch] = std::numeric_limits<double>::quiet_NaN();
        }
        
        // Fill MIP histograms
        if (!std::isnan(mMIPValues[ch])) {
            mHistMIPValues->SetBinContent(ch + 1, mMIPValues[ch]);
        }
        if (!std::isnan(mMIPDeviations[ch])) {
            mHistMIPDeviations->SetBinContent(ch + 1, mMIPDeviations[ch]);
        }
    }
    
    // Create MIP summary graph if it doesn't exist
    if (!mGraphMIPvsChannel) {
        mGraphMIPvsChannel = helper::registerGraph<TGraphErrors>(
            getObjectsManager(),
            quality_control::core::PublicationPolicy::ThroughStop,
            "AP", "MIPAnalysis/MIPVsChannel",
            "FT0: MIP values per channel;Channel ID;MIP value (ADC ch)",
            sNCHANNELS_PM);
            
        mGraphMIPvsChannel->GetXaxis()->SetLimits(-1, sNCHANNELS_PM + 1);
        mGraphMIPvsChannel->SetMarkerStyle(22);
        mGraphMIPvsChannel->SetMarkerSize(1.0);
        mGraphMIPvsChannel->SetMarkerColor(kBlue);
        mGraphMIPvsChannel->SetLineColor(kBlack);
        
        // Add expected MIP reference line
        auto* mipRefLine = new TLine(-0.5, mMIPExpectedValue, sNCHANNELS_PM + 0.5, mMIPExpectedValue);
        mipRefLine->SetLineColor(kGreen);
        mipRefLine->SetLineStyle(2);
        mipRefLine->SetLineWidth(2);
        mGraphMIPvsChannel->GetListOfFunctions()->Add(mipRefLine);
        
        ILOG(Info, Support) << "Created MIP analysis summary graph" << ENDM;
    }
    
    // Update MIP graph data
    for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
        mGraphMIPvsChannel->SetPoint(i, mChanX[i], mMIPValues[i]);
        mGraphMIPvsChannel->SetPointError(i, mChanXErr[i], 0.);
    }
}

void PostProcTask::updateMIPTrends()
{
    if (!mEnableMIPTracking || !mHistMIPTrends) return;
    
    // Get current time bin (simplified - use last bin)
    int timeBin = mHistMIPTrends->GetNbinsX();
    
    // Fill individual channel trends
    for (unsigned int ch = 0; ch < sNCHANNELS_PM; ++ch) {
        if (!std::isnan(mMIPDeviations[ch]) && ch < static_cast<unsigned int>(mHistMIPTrends->GetNbinsY())) {
            mHistMIPTrends->SetBinContent(timeBin, ch + 1, mMIPDeviations[ch]);
        }
    }
}

void PostProcTask::logAmplitudeStatistics() const
{
    if (mEnableAmplitudeAnalysis) {
        int totalFits = mSuccessfulFits + mFailedFits;
        if (totalFits > 0) {
            ILOG(Info, Support) << "Amplitude analysis: " << mSuccessfulFits << " successful fits ("
                               << (100.0 * mSuccessfulFits / totalFits) << "%), "
                               << mFailedFits << " failed fits ("
                               << (100.0 * mFailedFits / totalFits) << "%)" << ENDM;
        }
    }
    
    if (mEnableMIPTracking) {
        int totalMIPChannels = mMIPChannelsInRange + mMIPChannelsOutOfRange;
        if (totalMIPChannels > 0) {
            ILOG(Info, Support) << "MIP tracking: " << mMIPChannelsInRange << " channels in range ("
                               << (100.0 * mMIPChannelsInRange / totalMIPChannels) << "%), "
                               << mMIPChannelsOutOfRange << " out of range ("
                               << (100.0 * mMIPChannelsOutOfRange / totalMIPChannels) << "%)" << ENDM;
        }
    }
}

void PostProcTask::update(Trigger trg, framework::ServiceRegistryRef serviceReg)
{
  mPostProcHelper.update(trg, serviceReg);

  auto getBinContent2Ddiag = [](TH2F* hist, const std::string& binName) {
    const auto xBin = hist->GetXaxis()->FindBin(binName.c_str());
    const auto yBin = hist->GetYaxis()->FindBin(binName.c_str());
    return hist->GetBinContent(xBin, yBin);
  };

  // Trigger correlation
  auto hTrgCorr = mPostProcHelper.template getObject<TH2F>("TriggersCorrelation");
  mHistTriggers->Reset();

  if (hTrgCorr) {
    double totalStat{ 0 };
    for (int iBin = 1; iBin < mHistTriggers->GetXaxis()->GetNbins() + 1; iBin++) {
      const std::string binName{ mHistTriggers->GetXaxis()->GetBinLabel(iBin) };
      const auto binContent = getBinContent2Ddiag(hTrgCorr.get(), binName);
      mHistTriggers->SetBinContent(iBin, binContent);
      totalStat += binContent;
    }
    mHistTriggers->SetEntries(totalStat);
  }

  auto getVrtTrgCounters = [binPos = static_cast<int>(o2::ft0::Triggers::bitVertex) + 1](const auto& histTrgCnts) {
    return histTrgCnts->GetBinContent(binPos);
  };
  const auto trgVrtCnts = getVrtTrgCounters(mHistTriggers);
  // Trigger rates
  mHistTriggerRates->Reset();
  if (mPostProcHelper.IsNonEmptySample()) {
    mHistTriggerRates->Add(mHistTriggers.get());
    constexpr double factor = 1e3;                                                    // Hz -> kHz
    const double samplePeriod = 1. / (factor * mPostProcHelper.mCurrSampleLengthSec); // in sec^-1 units
    mHistTriggerRates->Scale(samplePeriod);
  }
  // PM bits
  auto hChDataBits = mPostProcHelper.template getObject<TH2F>("ChannelDataBits");
  auto hStatChannelID = mPostProcHelper.template getObject<TH1F>("StatChannelID");
  mHistChDataNOTbits->Reset();
  if (hChDataBits && hStatChannelID) {
    double totalStat{ 0 };
    for (int iBinX = 1; iBinX < hChDataBits->GetXaxis()->GetNbins() + 1; iBinX++) {
      for (int iBinY = 1; iBinY < hChDataBits->GetYaxis()->GetNbins() + 1; iBinY++) {
        const double nStatTotal = hStatChannelID->GetBinContent(iBinX);
        const double nStatPMbit = hChDataBits->GetBinContent(iBinX, iBinY);
        const double nStatNegPMbit = nStatTotal - nStatPMbit;
        totalStat += nStatNegPMbit;
        mHistChDataNOTbits->SetBinContent(iBinX, iBinY, nStatNegPMbit);
      }
    }
    mHistChDataNOTbits->SetEntries(totalStat);
  }
  // Amplitudes
  auto hAmpPerChannel = mPostProcHelper.template getObject<TH2F>("AmpPerChannel");
  if (hAmpPerChannel) {
    std::unique_ptr<TH1D> projNom(hAmpPerChannel->ProjectionX("projNom", hAmpPerChannel->GetYaxis()->FindBin(1.0), -1));
    std::unique_ptr<TH1D> projDen(hAmpPerChannel->ProjectionX("projDen"));
    mHistCFDEff->Divide(projNom.get(), projDen.get());

    // Reset amplitude histograms
    mHistAmpAInner->Reset();
    mHistAmpAOuter->Reset();
    mHistAmpC->Reset();
    mHistAmpNormPerChannel->Reset();
    if (mEnableMIPTracking) {
        mHistMIPValues->Reset();
        mHistMIPDeviations->Reset();
    }

    // Create projections for detector regions
    std::unique_ptr<TH1D> projAInner(hAmpPerChannel->ProjectionY("projAInner", 1, 32));
    std::unique_ptr<TH1D> projAOuter(hAmpPerChannel->ProjectionY("projAOuter", 33, 96));
    std::unique_ptr<TH1D> projC(hAmpPerChannel->ProjectionY("projC", 97, 208));

    mHistAmpAInner->Add(projAInner.get());
    mHistAmpAOuter->Add(projAOuter.get());
    mHistAmpC->Add(projC.get());

    // Normalize per channel
    int entries = mHistAmpNormPerChannel->GetEntries();
    for (int iBin = 1; iBin < hAmpPerChannel->GetXaxis()->GetNbins() + 1; iBin++) {
      std::unique_ptr<TH1D> ampForChannel(hAmpPerChannel->ProjectionY(Form("ampForChannel%i", iBin - 1), iBin, iBin));
      ampForChannel->Sumw2(kFALSE);

      const double integral = ampForChannel->Integral();
      const double scale = integral > 0. ? 1. / integral : 0.;
      mHistAmpNormPerChannel->Add(ampForChannel.get(), scale);
      entries += ampForChannel->GetEntries();
    }
    mHistAmpNormPerChannel->SetEntries(entries);
    
    // NEW: Perform enhanced amplitude analysis if enabled
    performAmplitudeAnalysis(hAmpPerChannel.get());
    
    // NEW: Perform MIP analysis if enabled
    if (mEnableAmplitudeAnalysis && mEnableMIPTracking) {
        performMIPAnalysis();
        updateMIPTrends();
    }
  }
  // Times
  auto hTimePerChannel = mPostProcHelper.template getObject<TH2F>("TimePerChannel");
  if (hTimePerChannel) {
    std::unique_ptr<TH1D> projInWindow(hTimePerChannel->ProjectionX("projInWindow", hTimePerChannel->GetYaxis()->FindBin(mLowTimeThreshold), hTimePerChannel->GetYaxis()->FindBin(mUpTimeThreshold)));
    std::unique_ptr<TH1D> projFull(hTimePerChannel->ProjectionX("projFull"));
    mHistTimeInWindow->Divide(projInWindow.get(), projFull.get());
    // hist with channel variables -> geometrical lot with channel variables
    const double scaleVrtTrg = trgVrtCnts > 0. ? 1. / trgVrtCnts : 0.0;
    projInWindow->Scale(scaleVrtTrg);
    mHistStatsSideA->Reset("content");
    mHistStatsSideC->Reset("content");
    mHistStatsSideA->SetStats(0);
    mHistStatsSideC->SetStats(0);
    mChannelGeometry.convertHist1D(projInWindow.get(), mHistStatsSideA.get(), mHistStatsSideC.get());
  }

  if (hAmpPerChannel && hTimePerChannel) {
    hAmpPerChannel->ProfileX("MeanAmplPerChannel");
    hTimePerChannel->ProfileX("MeanTimePerChannel");
    mAmpl->SetErrorOption("s");
    mTime->SetErrorOption("s");
    // for some reason the styling is not preserved after assigning result of ProfileX/Y() to already existing object
    mAmpl->SetMarkerStyle(8);
    mTime->SetMarkerStyle(8);
    mAmpl->SetLineColor(kBlack);
    mTime->SetLineColor(kBlack);
    mAmpl->SetDrawOption("P");
    mTime->SetDrawOption("P");
    mAmpl->GetXaxis()->SetTitleOffset(1);
    mTime->GetXaxis()->SetTitleOffset(1);
    mAmpl->GetYaxis()->SetTitleOffset(1);
    mTime->GetYaxis()->SetTitleOffset(1);
  }

  auto hBcVsTrg = mPostProcHelper.template getObject<TH2F>("BCvsTriggers");
  const auto& bcPattern = mPostProcHelper.getGRPLHCIFData().getBunchFilling();

  // Should be removed
  // BC pattern
  mHistBcPattern->Reset();
  for (int i = 0; i < sBCperOrbit + 1; i++) {
    for (int j = 0; j < mMapTechTrgBitsExtra.size() + 1; j++) {
      mHistBcPattern->SetBinContent(i + 1, j + 1, bcPattern.testBC(i));
    }
  }
  // Should be removed
  // Triggers in non-collision BCs
  mHistBcTrgOutOfBunchColl->Reset();
  float vmax = hBcVsTrg->GetBinContent(hBcVsTrg->GetMaximumBin());
  mHistBcTrgOutOfBunchColl->Add(hBcVsTrg.get(), mHistBcPattern.get(), 1, -1 * vmax);
  for (int i = 0; i < sBCperOrbit + 1; i++) {
    for (int j = 0; j < mMapTechTrgBitsExtra.size() + 1; j++) {
      if (mHistBcTrgOutOfBunchColl->GetBinContent(i + 1, j + 1) < 0) {
        mHistBcTrgOutOfBunchColl->SetBinContent(i + 1, j + 1, 0); // is it too slow?
      }
    }
  }
  mHistBcTrgOutOfBunchColl->SetEntries(mHistBcTrgOutOfBunchColl->Integral(1, sBCperOrbit, 1, mMapTechTrgBitsExtra.size()));
  for (int iBin = 1; iBin < mMapTechTrgBitsExtra.size() + 1; iBin++) {
    const std::string metadataKey = "BcVsTrgIntegralBin" + std::to_string(iBin);
    const std::string metadataValue = std::to_string(hBcVsTrg->Integral(1, sBCperOrbit, iBin, iBin));
    getObjectsManager()->getMonitorObject(mHistBcTrgOutOfBunchColl->GetName())->addOrUpdateMetadata(metadataKey, metadataValue);
    ILOG(Info, Support) << metadataKey << ":" << metadataValue << ENDM;
  }
  // Functor for in/out colliding BC fraction calculation
  auto calcFraction = [&pattern = bcPattern](const auto& histSrc, auto& histDst) {
    for (int iBin = 0; iBin < histSrc->GetYaxis()->GetNbins(); iBin++) {
      double cntInBC{ 0 };
      double cntOutOfBC{ 0 };
      const auto binPos = iBin + 1;
      std::unique_ptr<TH1D> proj(histSrc->ProjectionX("proj", binPos, binPos));
      for (int iBC = 0; iBC < proj->GetXaxis()->GetNbins(); iBC++) {
        if (pattern.testBC(iBC)) {
          cntInBC += proj->GetBinContent(iBC + 1);
        } else {
          cntOutOfBC += proj->GetBinContent(iBC + 1);
        }
      }
      const auto val = cntInBC > 0 ? cntOutOfBC / cntInBC : 0;
      histDst->SetBinContent(binPos, val);
    }
  };

  // New version for trigger fraction in non colliding BCa
  mHistTrg_outOfBC->Reset();
  calcFraction(hBcVsTrg, mHistTrg_outOfBC);

  // Channel stats in non-collision BCs
  auto hChIDvsBC = mPostProcHelper.template getObject<TH2F>("ChannelIDperBC");
  if (hChIDvsBC) {
    auto hChID_InBC = std::make_unique<TH1D>("hChID_InBC", "hChID_InBC", hChIDvsBC->GetYaxis()->GetNbins(), 0, hChIDvsBC->GetYaxis()->GetNbins());
    auto hChID_OutOfBC = std::make_unique<TH1D>("hChID_OutOfBC", "hChID_OutOfBC", hChIDvsBC->GetYaxis()->GetNbins(), 0, hChIDvsBC->GetYaxis()->GetNbins());
    if (mAsynchChannelLogic == std::string{ "standard" }) {
      // Standard logic, calculates outCollBC/inCollBC ratio
      calcFraction(hChIDvsBC, mHistChannelID_outOfBC);
    } else if (mAsynchChannelLogic == std::string{ "normalizedTrains" }) {
      // Train normlization, normalizes events in trains and then calculates outCollBC/inCollBC ratio
      const auto mapTrainBC = helper::getMapBCtrains(bcPattern.getBCPattern());
      for (int iChID = 0; iChID < hChIDvsBC->GetYaxis()->GetNbins(); iChID++) {
        double cntOutOfBC{ 0 };
        double cntInBC{ 0 };
        for (const auto& train : mapTrainBC) {
          double eventsTrain{ 0 };
          for (int iBC = train.first; iBC < train.first + train.second; iBC++) {
            eventsTrain += hChIDvsBC->GetBinContent(iBC + 1, iChID + 1);
          }
          if (train.second > 0) {
            cntInBC += eventsTrain / train.second;
          }
        }
        for (int iBC = 0; iBC < hChIDvsBC->GetXaxis()->GetNbins(); iBC++) {
          if (!bcPattern.testBC(iBC)) {
            cntOutOfBC += hChIDvsBC->GetBinContent(iBC + 1, iChID + 1);
          }
        }
        const auto val = cntInBC > 0 ? cntOutOfBC / cntInBC : 0;
        mHistChannelID_outOfBC->SetBinContent(iChID + 1, val);
      }
    }
  }

  // Fraction for trigger validation
  auto hTriggersSoftwareVsTCM = mPostProcHelper.template getObject<TH2F>("TriggersSoftwareVsTCM");
  if (hTriggersSoftwareVsTCM) {
    std::unique_ptr<TH1D> projOnlyHWorSW(hTriggersSoftwareVsTCM->ProjectionX("projOnlyHWorSW", 2, 3));
    std::unique_ptr<TH1D> projValidatedSWandHW(hTriggersSoftwareVsTCM->ProjectionX("projValidatedSWandHW", 4, 4));
    projOnlyHWorSW->LabelsDeflate();
    projValidatedSWandHW->LabelsDeflate();
    mHistTrgValidation->Divide(projOnlyHWorSW.get(), projValidatedSWandHW.get());
  }
  decomposeHists();
  logAmplitudeStatistics();
  setTimestampToMOs();
  // Needed for first-iter init
  mIsFirstIter = false;
  ILOG(Debug, Support) << "PostProcTask update completed with amplitude analysis" << ENDM;
}

void PostProcTask::decomposeHists()
{
  for (const auto& histName : mVecHistsToDecompose) {
    auto insertedMap = mMapHistsToDecompose.insert({ histName, {} });
    auto& mapHists = insertedMap.first->second;

    auto histSrcPtr = mPostProcHelper.template getObject<TH2F>(histName);
    if (histSrcPtr == nullptr) {
      continue;
    }
    const auto bins = histSrcPtr->GetYaxis()->GetNbins();
    const auto binLow = histSrcPtr->GetYaxis()->GetXmin();
    const auto binUp = histSrcPtr->GetYaxis()->GetXmax();

    for (const auto& chID : mVecChannelIDs) {
      auto insertedHistDst = mapHists.insert({ chID, nullptr });
      auto& histDstPtr = insertedHistDst.first->second;
      auto isInserted = insertedHistDst.second;
      if (isInserted == true) {
        // creation in first iter
        const std::string suffix = std::string{ Form("%03i", chID) };
        const std::string newHistName = histName + std::string{ "_" } + suffix;
        const std::string newHistTitle = histSrcPtr->GetTitle() + std::string{ " " } + suffix;
        histDstPtr = std::make_shared<HistDecomposed_t>(newHistName.c_str(), newHistTitle.c_str(), bins, binLow, binUp);
        getObjectsManager()->startPublishing(histDstPtr.get());
      }
      histDstPtr->Reset();
      // making projection
      const auto binPos = chID + 1;
      const std::unique_ptr<TH1D> proj(histSrcPtr->ProjectionY("proj", binPos, binPos));
      histDstPtr->Add(proj.get());
    }
  }
}

void PostProcTask::setTimestampToMOs()
{
  for (int iObj = 0; iObj < getObjectsManager()->getNumberPublishedObjects(); iObj++) {
    auto mo = getObjectsManager()->getMonitorObject(iObj);
    mo->addOrUpdateMetadata(mPostProcHelper.mTimestampMetaField, std::to_string(mPostProcHelper.mTimestampAnchor));
  }
}

void PostProcTask::finalize(Trigger t, framework::ServiceRegistryRef)
{
}

} // namespace o2::quality_control_modules::ft0
