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

// O2 QC / framework
#include "FT0/AmplitudePostProcTask.h"
#include "QualityControl/QcInfoLogger.h"

#include "FITCommon/HelperHist.h"
#include "FITCommon/HelperGraph.h"
#include "FITCommon/HelperCommon.h"

// ROOT
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <TList.h>
#include <TGraph.h>

// STL
#include <limits>
#include <memory>
#include <sstream>

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::fit;

namespace o2::quality_control_modules::ft0
{

// ---------------- configuration ----------------

void AmplitudePostProcTask::configure(const boost::property_tree::ptree& cfg)
{
  const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
  const char* cfgCustom = Form("%s.custom", cfgPath);

  mPostProcHelper.configure(cfg, cfgPath, "FT0");

  auto key = [&cfgCustom](const std::string& e) { return Form("%s.%s", cfgCustom, e.c_str()); };

  // ADC histogram config
  mAmpMin = helper::getConfigFromPropertyTree<int>(cfg, key("ampMin"), -100);
  mAmpMax = helper::getConfigFromPropertyTree<int>(cfg, key("ampMax"), 4100);
  mAmpBins = helper::getConfigFromPropertyTree<int>(cfg, key("ampBins"), 4200);

  // fixed fractional window
  mLeftSliceFrac = helper::getConfigFromPropertyTree<double>(cfg, key("leftSliceFrac"), 0.15);
  mRightSliceFrac = helper::getConfigFromPropertyTree<double>(cfg, key("rightSliceFrac"), 0.15);

  // weighted-mean bin cut
  mMinBinEntriesForWeight = helper::getConfigFromPropertyTree<int>(cfg, key("minBinEntriesForWeight"), 1);

  // expected gain (ADC/MIP)
  mExpectedGain = helper::getConfigFromPropertyTree<double>(cfg, key("expectedGain"), 14.0);

  // Digit task path
  mPathDigitQcTask = helper::getConfigFromPropertyTree<std::string>(cfg, key("pathDigitQcTask"), "FT0/MO/Digits/");

  // trend persistence knobs
  mTrendEnabled = helper::getConfigFromPropertyTree<bool>(cfg, key("trendEnabled"), true);
  mTrendScalarsFolder = helper::getConfigFromPropertyTree<std::string>(cfg, key("trendScalarsFolder"), "TrendsScalars");

  ILOG(Info, Support) << "AmplitudePostProcTask configured: "
                      << "amp=[" << mAmpMin << "," << mAmpMax << "] bins=" << mAmpBins
                      << ", window L/R=" << mLeftSliceFrac << "/" << mRightSliceFrac
                      << ", minBinEntriesForWeight=" << mMinBinEntriesForWeight
                      << ", expectedGain=" << mExpectedGain
                      << ", trendEnabled=" << (mTrendEnabled ? "yes" : "no")
                      << ", pathDigitQcTask=" << mPathDigitQcTask << ENDM;
}

// ---------------- lifecycle helpers ----------------

void AmplitudePostProcTask::reset()
{
  mMapHistAmpPerChannel.clear();

  mHistAmpAll.reset();
  mHistAmpAInner.reset();
  mHistAmpAOuter.reset();
  mHistAmpC.reset();
  mHistAmpNormPerChannel.reset();

  mGraphMeanNormVsChannel.reset();

  std::fill(mMeanW.begin(), mMeanW.end(), std::numeric_limits<double>::quiet_NaN());
  std::fill(mChanX.begin(), mChanX.end(), 0.);
  std::fill(mChanXErr.begin(), mChanXErr.end(), 0.);

  mMuAInner = mMuAOuter = mMuC = mMuAll = std::numeric_limits<double>::quiet_NaN();
  mSigAInner = mSigAOuter = mSigC = mSigAll = 0.;
}

void AmplitudePostProcTask::setTimestampToMOs()
{
  for (int iObj = 0; iObj < getObjectsManager()->getNumberPublishedObjects(); iObj++) {
    auto mo = getObjectsManager()->getMonitorObject(iObj);
    mo->addOrUpdateMetadata(mPostProcHelper.mTimestampMetaField,
                            std::to_string(mPostProcHelper.mTimestampAnchor));
  }
}

// ---------------- math & fitting ----------------

std::pair<double, double> AmplitudePostProcTask::computeWindow(double peak) const
{
  // symmetric or asymmetric fractional window around the peak
  double xmin = std::max<double>(peak - mLeftSliceFrac * std::abs(peak), mAmpMin);
  double xmax = std::min<double>(peak + mRightSliceFrac * std::abs(peak), mAmpMax);
  if (xmax <= xmin) {
    // fallback to tiny window around peak
    xmin = std::max<double>(peak - 1.0, mAmpMin);
    xmax = std::min<double>(peak + 1.0, mAmpMax);
  }
  return { xmin, xmax };
}

double AmplitudePostProcTask::computeWeightedMeanInWindow(const TH1D* h, double xmin, double xmax) const
{
  if (!h) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const int bmin = std::max(1, h->GetXaxis()->FindBin(xmin));
  const int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(xmax));

  double wsum = 0.;
  double xwsum = 0.;

  for (int b = bmin; b <= bmax; ++b) {
    const double n = h->GetBinContent(b);
    if (n < mMinBinEntriesForWeight)
      continue; // ignore empty/tiny bins
    const double x = h->GetBinCenter(b);
    wsum += n;
    xwsum += x * n;
  }

  if (wsum <= 0.)
    return std::numeric_limits<double>::quiet_NaN();
  return xwsum / wsum;
}

bool AmplitudePostProcTask::fitRegionGaussianAndOverlay(TH1F* regionHist,
                                                        const char* fitName,
                                                        double& outMu,
                                                        double& outSigma) const
{
  if (!regionHist || regionHist->GetEntries() < 50) {
    outMu = std::numeric_limits<double>::quiet_NaN();
    outSigma = 0.;
    return false;
  }

  // 1) Find peak and the window around it
  const int bMax = regionHist->GetMaximumBin();
  const double peak = regionHist->GetBinCenter(bMax);
  const auto [xmin, xmax] = computeWindow(peak);

  // 2) SAFELY remove any previous overlays by name - let ROOT manage memory
  auto* lst = regionHist->GetListOfFunctions();

  // 3) Fit a Gaussian, but DO NOT attach TF1 to the histogram
  TF1 fG_local("fG_tmp", "gaus", xmin, xmax);
  regionHist->Fit(&fG_local, "QNR", "", xmin, xmax); // Q: quiet, N: don't store, R: range

  outMu = fG_local.GetParameter(1);
  outSigma = std::abs(fG_local.GetParameter(2));

  // 4) Create a TGraph curve from the fit parameters and attach that instead of TF1
  {
    const int n = 200;
    std::vector<double> xs(n), ys(n);
    const double dx = (xmax - xmin) / (n - 1);
    for (int i = 0; i < n; ++i) {
      double x = xmin + i * dx;
      xs[i] = x;
      ys[i] = fG_local.Eval(x);
    }
    auto* gr = new TGraph(n, xs.data(), ys.data());
    gr->SetName(Form("%s_curve", fitName));
    gr->SetTitle("Gaussian fit");
    gr->SetLineColor(kGreen + 2);
    gr->SetLineWidth(2);
    if (lst)
      lst->Add(gr); // ROOT takes ownership
  }

  // 5) Weighted mean in the same window (red vertical line)
  double muW = std::numeric_limits<double>::quiet_NaN();
  {
    double wsum = 0., xwsum = 0.;
    const int bmin = std::max(1, regionHist->GetXaxis()->FindBin(xmin));
    const int bmax = std::min(regionHist->GetNbinsX(), regionHist->GetXaxis()->FindBin(xmax));
    for (int b = bmin; b <= bmax; ++b) {
      const double n = regionHist->GetBinContent(b);
      if (n < mMinBinEntriesForWeight)
        continue;
      const double x = regionHist->GetBinCenter(b);
      wsum += n;
      xwsum += x * n;
    }
    if (wsum > 0.)
      muW = xwsum / wsum;
  }
  if (!std::isnan(muW)) {
    auto* wm = new TLine(muW, 0., muW, regionHist->GetMaximum());
    wm->SetLineColor(kRed);
    wm->SetLineStyle(1);
    wm->SetLineWidth(2);
    if (lst)
      lst->Add(wm); // ROOT takes ownership
  }

  // 6) Window bounds (blue dashed lines)
  auto* lwmin = new TLine(xmin, 0., xmin, regionHist->GetMaximum());
  lwmin->SetLineColor(kBlue);
  lwmin->SetLineStyle(2);
  lwmin->SetLineWidth(2);
  if (lst)
    lst->Add(lwmin); // ROOT takes ownership

  auto* lwmax = new TLine(xmax, 0., xmax, regionHist->GetMaximum());
  lwmax->SetLineColor(kBlue);
  lwmax->SetLineStyle(2);
  lwmax->SetLineWidth(2);
  if (lst)
    lst->Add(lwmax); // ROOT takes ownership

  // 7) Basic sanity
  if (outMu < mAmpMin || outMu > mAmpMax || outSigma <= 0.) {
    ILOG(Warning, Support) << "Suspicious region fit: mu=" << outMu << " sigma=" << outSigma << ENDM;
    return false;
  }
  return true;
}

// ---------------- initialize ----------------

void AmplitudePostProcTask::initialize(Trigger trig, framework::ServiceRegistryRef services)
{
  ILOG(Info, Support) << "Initializing AmplitudePostProcTask (FT0)" << ENDM;

  mPostProcHelper.initialize(trig, services);
  reset();

  // global histos
  mHistAmpAll = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpAllChannels", "FT0 all channels;Channel amplitude (ADC ch);Counts",
    mAmpBins, mAmpMin, mAmpMax);

  mHistAmpAInner = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpAInner", "FT0 A-side inner channels (0-31);Channel amplitude (ADC ch);Counts",
    mAmpBins, mAmpMin, mAmpMax);

  mHistAmpAOuter = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpAOuter", "FT0 A-side outer channels (32-95);Channel amplitude (ADC ch);Counts",
    mAmpBins, mAmpMin, mAmpMax);

  mHistAmpC = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpC", "FT0 C-side channels (96-207);Channel amplitude (ADC ch);Counts",
    mAmpBins, mAmpMin, mAmpMax);

  mHistAmpNormPerChannel = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "AmpNormPerChannel", "FT0 normalized per channel;Channel amplitude (ADC ch);Normalized counts",
    mAmpBins, mAmpMin, mAmpMax);
  mHistAmpNormPerChannel->Sumw2(kFALSE);

  mTrendAInner = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "TrendsScalars/AInner", "AInner;dummy;value", 1, mAmpMin, mAmpMax);
  mTrendAOuter = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "TrendsScalars/AOuter", "AOuter;dummy;value", 1, mAmpMin, mAmpMax);
  mTrendC = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "TrendsScalars/C", "C;dummy;value", 1, mAmpMin, mAmpMax);
  mTrendAll = helper::registerHist<TH1F>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "", "TrendsScalars/All", "All;dummy;value", 1, mAmpMin, mAmpMax);

  // per-channel histos (0..207)
  for (unsigned int ch = 0; ch < sNCHANNELS_PM; ++ch) {
    const std::string name = Form("AmplitudePerChannel/Amp_ch%03u", ch);
    const std::string title = Form("FT0 channel %u;Channel amplitude (ADC ch);Counts", ch);
    mMapHistAmpPerChannel[ch] = helper::registerHist<TH1F>(
      getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
      "", name, title,
      mAmpBins, mAmpMin, mAmpMax);
  }

  // per-channel summary graph: muW / expectedGain vs channel
  mGraphMeanNormVsChannel = helper::registerGraph<TGraphErrors>(
    getObjectsManager(), quality_control::core::PublicationPolicy::ThroughStop,
    "AP", "GaussianSummary/MeanNormVsChannel",
    "FT0: Weighted mean normalized;Channel ID;#mu_{W}/expectedGain", sNCHANNELS_PM);
  mGraphMeanNormVsChannel->GetXaxis()->SetLimits(-1, sNCHANNELS_PM + 1);
  mGraphMeanNormVsChannel->GetYaxis()->SetRangeUser(0., 2.);
  mGraphMeanNormVsChannel->SetMarkerStyle(21);
  mGraphMeanNormVsChannel->SetMarkerSize(1.1);
  mGraphMeanNormVsChannel->SetMarkerColor(kRed);
  mGraphMeanNormVsChannel->SetLineColor(kBlack);
  {
    auto* ref = new TLine(-0.5, 1.0, sNCHANNELS_PM + 0.5, 1.0);
    ref->SetLineColor(kBlue);
    ref->SetLineStyle(2);
    ref->SetLineWidth(2);
    mGraphMeanNormVsChannel->GetListOfFunctions()->Clear();
    mGraphMeanNormVsChannel->GetListOfFunctions()->Add(ref);
  }

  ILOG(Info, Support) << "Initialization complete" << ENDM;
}

// ---------------- update ----------------

void AmplitudePostProcTask::update(Trigger trig, framework::ServiceRegistryRef serviceReg)
{
  ILOG(Info, Support) << "=== AmplitudePostProcTask::update() START ===" << ENDM;

  try {
    mPostProcHelper.update(trig, serviceReg);
    ILOG(Info, Support) << "PostProcHelper updated, timestamp = " << mPostProcHelper.mTimestampAnchor << ENDM;

    // fetch 2D amplitude per channel from the Digit task
    auto h2 = mPostProcHelper.template getObject<TH2F>("AmpPerChannel");
    if (!h2) {
      ILOG(Error, Support) << "MO 'AmpPerChannel' not found at path " << mPathDigitQcTask << ENDM;
      setTimestampToMOs();
      return;
    }
    ILOG(Info, Support) << "Successfully retrieved AmpPerChannel histogram" << ENDM;

    // reset aggregates
    mHistAmpAll->Reset();
    mHistAmpAInner->Reset();
    mHistAmpAOuter->Reset();
    mHistAmpC->Reset();
    mHistAmpNormPerChannel->Reset();

    // clear region fit results
    mMuAInner = mMuAOuter = mMuC = mMuAll = std::numeric_limits<double>::quiet_NaN();
    mSigAInner = mSigAOuter = mSigC = mSigAll = 0.;

    // iterate channels
    for (int chBin = 1, nBins = h2->GetXaxis()->GetNbins(); chBin <= nBins; ++chBin) {
      unsigned int ch = chBin - 1;
      if (ch >= sNCHANNELS_PM)
        continue;

      std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("p_ch%03u", ch), chBin, chBin));
      proj->Sumw2(kFALSE);

      // add to per-channel and region aggregates
      mMapHistAmpPerChannel[ch]->Reset();
      mMapHistAmpPerChannel[ch]->Add(proj.get());
      mHistAmpAll->Add(proj.get());

      if (ch < 32) {
        mHistAmpAInner->Add(proj.get());
      } else if (ch < 96) {
        mHistAmpAOuter->Add(proj.get());
      } else {
        mHistAmpC->Add(proj.get());
      }

      // normalized per-channel stack
      const double intg = proj->Integral();
      if (intg > 0.) {
        mHistAmpNormPerChannel->Add(proj.get(), 1.0 / intg);
      }

      // per-channel weighted mean in window
      if (proj->GetEntries() > 0) {
        const int bMax = proj->GetMaximumBin();
        const double peak = proj->GetBinCenter(bMax);
        const auto [xmin, xmax] = computeWindow(peak);
        mMeanW[ch] = computeWeightedMeanInWindow(proj.get(), xmin, xmax);
      } else {
        mMeanW[ch] = std::numeric_limits<double>::quiet_NaN();
      }

      mChanX[ch] = ch;
      mChanXErr[ch] = 0.;
    }

    // region gaussian fits + overlays
    fitRegionGaussianAndOverlay(mHistAmpAInner.get(), "fG_AInner", mMuAInner, mSigAInner);
    fitRegionGaussianAndOverlay(mHistAmpAOuter.get(), "fG_AOuter", mMuAOuter, mSigAOuter);
    fitRegionGaussianAndOverlay(mHistAmpC.get(), "fG_C", mMuC, mSigC);
    fitRegionGaussianAndOverlay(mHistAmpAll.get(), "fG_All", mMuAll, mSigAll);

    // update per-channel summary graph (muW / expectedGain)
    for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
      const double y = (std::isnan(mMeanW[i]) || mExpectedGain <= 0.)
                         ? std::numeric_limits<double>::quiet_NaN()
                         : (mMeanW[i] / mExpectedGain);
      mGraphMeanNormVsChannel->SetPoint(i, mChanX[i], y);
      mGraphMeanNormVsChannel->SetPointError(i, mChanXErr[i], 0.);
    }

    auto updateHistogram = [](TH1F* h, double v) {
      if (!h || std::isnan(v))
        return;
      h->Reset("ICES");
      h->Fill(v);
    };

    updateHistogram(mTrendAInner.get(), mMuAInner);
    updateHistogram(mTrendAOuter.get(), mMuAOuter);
    updateHistogram(mTrendC.get(), mMuC);
    updateHistogram(mTrendAll.get(), mMuAll);

    setTimestampToMOs();

    ILOG(Info, Support) << "Update done. Region mu (ADC): "
                        << "AInner=" << mMuAInner
                        << ", AOuter=" << mMuAOuter
                        << ", C=" << mMuC
                        << ", All=" << mMuAll << ENDM;

  } catch (const std::exception& e) {
    ILOG(Error, Support) << "Exception in AmplitudePostProcTask::update(): " << e.what() << ENDM;
    throw;
  } catch (...) {
    ILOG(Error, Support) << "Unknown exception in AmplitudePostProcTask::update()" << ENDM;
    throw;
  }

  ILOG(Info, Support) << "=== AmplitudePostProcTask::update() END ===" << ENDM;
}

// Finalize can be empty now
void AmplitudePostProcTask::finalize(Trigger, framework::ServiceRegistryRef)
{
  ILOG(Info, Support) << "Finalize - trend objects already exist and updated" << ENDM;
}

} // namespace o2::quality_control_modules::ft0