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
#include "QualityControl/MonitorObject.h"
#include "Common/Utils.h"
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
        mCcdbUrl = cfg.get_child("qc.config.conditionDB.url").get_value<std::string>();

        const char* cfgPath   = Form("qc.postprocessing.%s", getID().c_str());
        const char* cfgCustom = Form("%s.custom", cfgPath);
        auto full = [&cfgCustom](std::string_view entry) {
            return Form("%s.%s", cfgCustom, std::string(entry).c_str());
        };

        mPathDigitQcTask    = helper::getConfigFromPropertyTree<std::string>(
            cfg, full("pathDigitQcTask"), "FV0/MO/DigitsNew");
        mTimestampMetaField = helper::getConfigFromPropertyTree<std::string>(
            cfg, full("timestampMetaField"), "timestampTF");
        mAmpMin    = helper::getConfigFromPropertyTree<int >(cfg, full("ampMin"),   -100);
        mAmpMax    = helper::getConfigFromPropertyTree<int >(cfg, full("ampMax"),   4100);
        mAmpBins   = helper::getConfigFromPropertyTree<int >(cfg, full("ampBins"),  4200);
        mSliceFrac = helper::getConfigFromPropertyTree<double>(cfg, full("sliceFrac"),0.25);

        ILOG(Info, Support) << "AmplitudePostProcTask configured" << ENDM;
    }
    
    void AmplitudePostProcTask::reset()
    {
        mMapHistAmpPerChannel.clear();
        mHistAmpAll.reset();
        mHistAmpNormPerChannel.reset();
    }

    void AmplitudePostProcTask::setTimestampToMOs(long long ts)
    {
        for (int i = 0, n = getObjectsManager()->getNumberPublishedObjects();
            i < n; ++i) {
            if (auto* mo = getObjectsManager()->getMonitorObject(i)) {
                mo->addOrUpdateMetadata(mTimestampMetaField, std::to_string(ts));
            }
        }
    }

    void AmplitudePostProcTask::initialize(Trigger,
                                        framework::ServiceRegistryRef services)
    {
        ILOG(Info, Support) << "Initialising AmplitudePostProcTask" << ENDM;
        reset();

        mDatabase = &services.get<o2::quality_control::repository::DatabaseInterface>();
        mCcdbApi.init(mCcdbUrl);

        // Global histograms
        mHistAmpAll = helper::registerHist<TH1F>(getObjectsManager(),
        quality_control::core::PublicationPolicy::ThroughStop,
        "", "AmpAllChannels",
        "FV0 all channels;Channel amplitude (ADC ch);Counts",
        mAmpBins, mAmpMin, mAmpMax);

        mHistAmpNormPerChannel = helper::registerHist<TH1F>(getObjectsManager(),
        quality_control::core::PublicationPolicy::ThroughStop,
        "", "AmpNormPerChannel",
        "FV0 normalised (per channel);Channel amplitude (ADC ch);Normalised counts",
        mAmpBins, mAmpMin, mAmpMax);
        mHistAmpNormPerChannel->Sumw2(kFALSE);

        // Per-channel histograms
        for (unsigned int ch = 0; ch < sNCHANNELS_PM; ++ch) {
            const std::string name  = Form("AmplitudePerChannel/Amp_ch%02u", ch);
            const std::string title = Form("FV0 channel %u;Channel amplitude (ADC ch);Counts", ch);

            mMapHistAmpPerChannel[ch] = helper::registerHist<TH1F>(getObjectsManager(),
            quality_control::core::PublicationPolicy::ThroughStop,
            "", name, title,
            mAmpBins, mAmpMin, mAmpMax);
        }
    }

    void AmplitudePostProcTask::update(Trigger trig, framework::ServiceRegistryRef)
    {
        constexpr double kScale = 16.;  // pp scaling for now

        // Retrieve source histogram
        auto mo = mDatabase->retrieveMO(mPathDigitQcTask,
                                    "AmpPerChannel", trig.timestamp, trig.activity);
        auto* h2 = mo ? dynamic_cast<TH2F*>(mo->getObject()) : nullptr;
        if (!h2) {
            ILOG(Error, Support) << "MO 'AmpPerChannel' not found at " << mPathDigitQcTask << ENDM;
            return;
        }

        // Reset summary hists
        mHistAmpAll->Reset();
        mHistAmpNormPerChannel->Reset();
        for (auto& [_, h] : mMapHistAmpPerChannel) h->Reset();

        // Channel loop
        for (int chBin = 1, nBins = h2->GetXaxis()->GetNbins(); chBin <= nBins; ++chBin) {
            unsigned int ch = chBin - 1;
            if (ch >= sNCHANNELS_PM) continue;

            std::unique_ptr<TH1D> proj(h2->ProjectionY(Form("p_ch%02u", ch), chBin, chBin));
            proj->Sumw2(kFALSE);

            mMapHistAmpPerChannel[ch]->Add(proj.get());
            mHistAmpAll->Add(proj.get());

            if (double intg = proj->Integral(); intg > 0.)
                mHistAmpNormPerChannel->Add(proj.get(), 1. / intg);

            // Gaussian slice fit
            if (proj->GetEntries() > 50) {
                static TF1 fG("fG", "gaus", mAmpMin, mAmpMax);
                int    bMax = proj->GetMaximumBin();
                double peak = proj->GetBinCenter(bMax);
                double xmin = std::max<double>(peak - mSliceFrac * std::abs(peak), mAmpMin);
                double xmax = std::min<double>(peak + mSliceFrac * std::abs(peak), mAmpMax);

                proj->Fit(&fG, "QNR", "", xmin, xmax);

                mMean [ch] = fG.GetParameter(1);
                mSigma[ch] = std::abs(fG.GetParameter(2));
            } else {
                mMean [ch] = std::numeric_limits<double>::quiet_NaN();
                mSigma[ch] = 0.;
            }
            mChanX    [ch] = ch;
            mChanXErr [ch] = 0.;
        }

        // Build / refresh MPV graph
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
            mGraphMPVDiv16->SetMarkerSize (1.1);
            mGraphMPVDiv16->SetMarkerColor(kRed);
            mGraphMPVDiv16->SetLineColor  (kBlack);

            // Add reference line
            mGraphMPVDiv16->GetListOfFunctions()->Clear();
            auto* refLine = new TLine(-0.5, 1.0, 48.5, 1.0);
            refLine->SetLineColor(kBlue);
            refLine->SetLineStyle(2);
            refLine->SetLineWidth(2);
            mGraphMPVDiv16->GetListOfFunctions()->Add(refLine);
        }

        for (std::size_t i = 0; i < sNCHANNELS_PM; ++i) {
            mGraphMPVDiv16->SetPoint      (i, mChanX[i],    mMean[i]  / kScale);
            mGraphMPVDiv16->SetPointError (i, mChanXErr[i], mSigma[i] / kScale);
        }

        // Propagate time stamp
        long long ts = trig.timestamp;
        if (mo) {
            auto it = mo->getMetadataMap().find(mTimestampMetaField);
            if (it != mo->getMetadataMap().end()) {
                try { ts = std::stoll(it->second); } 
                catch (...) {
                    ILOG(Error, Support) << "Failed to parse timestamp from metadata" << ENDM;
                }
            }
        }
        setTimestampToMOs(ts);
    }

    void AmplitudePostProcTask::finalize(Trigger, framework::ServiceRegistryRef)
    {
        ILOG(Info, Support) << "AmplitudePostProcTask finalised" << ENDM;
    }
}   // namespace o2::quality_control_modules::fv0