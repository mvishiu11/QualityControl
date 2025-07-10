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

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::fit;

namespace o2::quality_control_modules::fv0
{

void AmplitudePostProcTask::configure(const boost::property_tree::ptree& cfg)
{
    // Configure PostProcHelper
    const char* cfgPath = Form("qc.postprocessing.%s", getID().c_str());
    const char* cfgCustom = Form("%s.custom", cfgPath);
    
    // Initialize PostProcHelper with FV0 detector name
    mPostProcHelper.configure(cfg, cfgPath, "FV0");
    
    // Configure amplitude-specific parameters
    auto cfgGet = [&cfgCustom](const std::string& entry) {
        return Form("%s.%s", cfgCustom, entry.c_str());
    };
    
    mAmpMin = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampMin"), -100);
    mAmpMax = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampMax"), 4100);
    mAmpBins = helper::getConfigFromPropertyTree<int>(cfg, cfgGet("ampBins"), 4200);
    mSliceFrac = helper::getConfigFromPropertyTree<double>(cfg, cfgGet("sliceFrac"), 0.25);
    
    ILOG(Info, Support) << "AmplitudePostProcTask configured with slice fraction: " << mSliceFrac << ENDM;
}

void AmplitudePostProcTask::reset()
{
    mMapHistAmpPerChannel.clear();
    mHistAmpAll.reset();
    mHistAmpNormPerChannel.reset();
    mGraphMPVDiv16.reset();
}

void AmplitudePostProcTask::setTimestampToMOs()
{
    for (int iObj = 0; iObj < getObjectsManager()->getNumberPublishedObjects(); iObj++) {
        auto mo = getObjectsManager()->getMonitorObject(iObj);
        mo->addOrUpdateMetadata(mPostProcHelper.mTimestampMetaField, 
                               std::to_string(mPostProcHelper.mTimestampAnchor));
    }
}

void AmplitudePostProcTask::initialize(Trigger trig, framework::ServiceRegistryRef services)
{
    ILOG(Info, Support) << "Initialising AmplitudePostProcTask" << ENDM;
    mPostProcHelper.initialize(trig, services);
    reset();
    
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
    
    ILOG(Info, Support) << "AmplitudePostProcTask initialized with " << sNCHANNELS_PM << " channels" << ENDM;
}

void AmplitudePostProcTask::update(Trigger trig, framework::ServiceRegistryRef serviceReg)
{
    mPostProcHelper.update(trig, serviceReg);
    
    constexpr double kScale = 16.0;  // pp scaling factor
    
    auto h2 = mPostProcHelper.template getObject<TH2F>("AmpPerChannel");
    if (!h2) {
        ILOG(Error, Support) << "MO 'AmpPerChannel' not found in PostProcHelper" << ENDM;
        setTimestampToMOs();
        return;
    }
    
    ILOG(Debug, Support) << "Processing amplitude data with " << h2->GetEntries() << " entries" << ENDM;
    
    mHistAmpAll->Reset();
    mHistAmpNormPerChannel->Reset();
    for (auto& [_, h] : mMapHistAmpPerChannel) {
        h->Reset();
    }
    
    // Process each channel
    for (int chBin = 1, nBins = h2->GetXaxis()->GetNbins(); chBin <= nBins; ++chBin) {
        unsigned int ch = chBin - 1;
        if (ch >= sNCHANNELS_PM) continue;
        
        std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("p_ch%02u", ch), chBin, chBin));
        proj->Sumw2(kFALSE);

        mMapHistAmpPerChannel[ch]->Add(proj.get());
        mHistAmpAll->Add(proj.get());

        if (double intg = proj->Integral(); intg > 0.) {
            mHistAmpNormPerChannel->Add(proj.get(), 1.0 / intg);
        }
        
        // Perform Gaussian slice fit for peak analysis
        if (proj->GetEntries() > 50) {
            static TF1 fG("fG", "gaus", mAmpMin, mAmpMax);
            
            int bMax = proj->GetMaximumBin();
            double peak = proj->GetBinCenter(bMax);
            
            double xmin = std::max<double>(peak - mSliceFrac * std::abs(peak), mAmpMin);
            double xmax = std::min<double>(peak + mSliceFrac * std::abs(peak), mAmpMax);

            proj->Fit(&fG, "QNR", "", xmin, xmax);

            mMean[ch] = fG.GetParameter(1);             // Gaussian mean
            mSigma[ch] = std::abs(fG.GetParameter(2));  // Gaussian sigma
        } else {
            // Insufficient statistics for reliable fit
            mMean[ch] = std::numeric_limits<double>::quiet_NaN();
            mSigma[ch] = 0.0;
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
            
        // Configure graph appearance
        mGraphMPVDiv16->GetXaxis()->SetLimits(-1, 49);
        mGraphMPVDiv16->GetYaxis()->SetRangeUser(0., 2.);
        mGraphMPVDiv16->SetMarkerStyle(21);
        mGraphMPVDiv16->SetMarkerSize(1.1);
        mGraphMPVDiv16->SetMarkerColor(kRed);
        mGraphMPVDiv16->SetLineColor(kBlack);
        
        // Add reference line at y=1.0
        mGraphMPVDiv16->GetListOfFunctions()->Clear();
        auto* refLine = new TLine(-0.5, 1.0, 48.5, 1.0);
        refLine->SetLineColor(kBlue);
        refLine->SetLineStyle(2);
        refLine->SetLineWidth(2);
        mGraphMPVDiv16->GetListOfFunctions()->Add(refLine);
        
        ILOG(Info, Support) << "Created Gaussian summary graph with reference line" << ENDM;
    }
    
    // Update graph points with scaled fit results
    for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
        mGraphMPVDiv16->SetPoint(i, mChanX[i], mMean[i] / kScale);
        mGraphMPVDiv16->SetPointError(i, mChanXErr[i], mSigma[i] / kScale);
    }
    
    setTimestampToMOs();
    
    ILOG(Info, Support) << "AmplitudePostProcTask update completed successfully" << ENDM;
}

void AmplitudePostProcTask::finalize(Trigger, framework::ServiceRegistryRef)
{
    int totalChannelsProcessed = 0;
    for (auto& [ch, hist] : mMapHistAmpPerChannel) {
        if (hist && hist->GetEntries() > 0) {
            totalChannelsProcessed++;
            ILOG(Debug, Support) << "Channel " << ch << " processed " 
                               << hist->GetEntries() << " entries" << ENDM;
        }
    }
    
    if (mGraphMPVDiv16) {
        ILOG(Info, Support) << "Gaussian summary graph contains " 
                           << mGraphMPVDiv16->GetN() << " points" << ENDM;
    }
    
    ILOG(Info, Support) << "AmplitudePostProcTask finalized. Processed " 
                       << totalChannelsProcessed << " channels." << ENDM;
}

} // namespace o2::quality_control_modules::fv0