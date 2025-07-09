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

#include "FV0/AmplitudePostProcTask.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/MonitorObject.h"
#include "Common/Utils.h"
#include "FITCommon/HelperHist.h"
#include "FITCommon/HelperCommon.h"

#include <TH1F.h>
#include <TH2F.h>
#include <memory>

using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::fit;

namespace o2::quality_control_modules::fv0
{

AmplitudePostProcTask::~AmplitudePostProcTask()
{
  reset();
}

void AmplitudePostProcTask::configure(const boost::property_tree::ptree& config)
{
  mCcdbUrl = config.get_child("qc.config.conditionDB.url").get_value<std::string>();
  
  const char* configPath = Form("qc.postprocessing.%s", getID().c_str());
  const char* configCustom = Form("%s.custom", configPath);
  
  ILOG(Info, Support) << "Configuring AmplitudePostProcTask with path: " << configPath << ENDM;
  
  auto cfgPath = [&configCustom](const std::string& entry) {
    return Form("%s.%s", configCustom, entry.c_str());
  };
  
  // Configure source task path
  auto node = config.get_child_optional(Form("%s.custom.pathDigitQcTask", configPath));
  if (node) {
    mPathDigitQcTask = node.get_ptr()->get_child("").get_value<std::string>();
    ILOG(Debug, Support) << "Using configured pathDigitQcTask: " << mPathDigitQcTask << ENDM;
  } else {
    mPathDigitQcTask = "FV0/MO/DigitsNew";
    ILOG(Debug, Support) << "Using default pathDigitQcTask: " << mPathDigitQcTask << ENDM;
  }
  
  // Configure histogram parameters
  mTimestampMetaField = helper::getConfigFromPropertyTree<std::string>(config, cfgPath("timestampMetaField"), "timestampTF");
  mAmpMin = helper::getConfigFromPropertyTree<int>(config, cfgPath("ampMin"), -100);
  mAmpMax = helper::getConfigFromPropertyTree<int>(config, cfgPath("ampMax"), 4100);
  mAmpBins = helper::getConfigFromPropertyTree<int>(config, cfgPath("ampBins"), 4200);
  
  ILOG(Info, Support) << "Amplitude histogram configuration: " << mAmpBins 
                      << " bins, range [" << mAmpMin << ", " << mAmpMax << "]" << ENDM;
}

void AmplitudePostProcTask::reset()
{
  mMapHistAmpPerChannel.clear();
  mHistAmpAll.reset();
  mHistAmpNormPerChannel.reset();
}

void AmplitudePostProcTask::initialize(Trigger, framework::ServiceRegistryRef services)
{
  ILOG(Info, Support) << "Initializing AmplitudePostProcTask" << ENDM;
  
  reset();
  
  // Initialize database interface
  mDatabase = &services.get<o2::quality_control::repository::DatabaseInterface>();
  mCcdbApi.init(mCcdbUrl);
  
  // Create combined amplitude histogram
  mHistAmpAll = helper::registerHist<TH1F>(
    getObjectsManager(), 
    quality_control::core::PublicationPolicy::ThroughStop, 
    "", 
    "AmpAllChannels", 
    "FV0 All Channel Amplitudes;Channel amplitude (ADC ch);Counts", 
    mAmpBins, mAmpMin, mAmpMax
  );
  
  // Create normalized amplitude histogram
  mHistAmpNormPerChannel = helper::registerHist<TH1F>(
    getObjectsManager(), 
    quality_control::core::PublicationPolicy::ThroughStop, 
    "", 
    "AmpNormPerChannel", 
    "FV0 Normalized Channel Amplitudes;Channel amplitude (ADC ch);Normalized counts", 
    mAmpBins, mAmpMin, mAmpMax
  );
  mHistAmpNormPerChannel->Sumw2(kFALSE);
  
  // Create individual histograms for each channel
  for (unsigned int chID = 0; chID < sNCHANNELS_PM; chID++) {
    const std::string histName = Form("AmplitudePerChannel/Amp_channel%03d", chID);
    const std::string histTitle = Form("FV0 Amplitude Channel %d;Channel amplitude (ADC ch);Counts", chID);
    
    mMapHistAmpPerChannel[chID] = helper::registerHist<TH1F>(
      getObjectsManager(), 
      quality_control::core::PublicationPolicy::ThroughStop, 
      "", 
      histName, 
      histTitle, 
      mAmpBins, mAmpMin, mAmpMax
    );
  }
  
  ILOG(Info, Support) << "Initialized " << sNCHANNELS_PM << " channel histograms" << ENDM;
}

void AmplitudePostProcTask::update(Trigger t, framework::ServiceRegistryRef)
{
  ILOG(Debug, Devel) << "Updating AmplitudePostProcTask" << ENDM;
  
  // Retrieve the 2D amplitude vs channel histogram from digit QC task
  auto mo = mDatabase->retrieveMO(mPathDigitQcTask, "AmpPerChannel", t.timestamp, t.activity);
  auto hAmpPerChannel = mo ? dynamic_cast<TH2F*>(mo->getObject()) : nullptr;
  
  if (!hAmpPerChannel) {
    ILOG(Error, Support) << "Failed to retrieve MO \"AmpPerChannel\" from " << mPathDigitQcTask << ENDM;
    return;
  }
  
  // Reset all histograms before filling
  mHistAmpAll->Reset();
  mHistAmpNormPerChannel->Reset();
  
  for (auto& [chID, hist] : mMapHistAmpPerChannel) {
    hist->Reset();
  }
  
  // Process each channel
  const int nChannelBins = hAmpPerChannel->GetXaxis()->GetNbins();
  ILOG(Debug, Support) << "Processing " << nChannelBins << " channel bins" << ENDM;
  
  for (int chBin = 1; chBin <= nChannelBins; chBin++) {
    const unsigned int chID = chBin - 1;
    
    if (chID >= sNCHANNELS_PM) {
      ILOG(Warning, Support) << "Channel bin " << chBin << " (chID=" << chID 
                            << ") exceeds maximum channels (" << sNCHANNELS_PM << ")" << ENDM;
      continue;
    }
    
    // Create projection for this channel
    std::unique_ptr<TH1D> projChannel(
      hAmpPerChannel->ProjectionY(Form("proj_ch%d", chID), chBin, chBin)
    );
    projChannel->Sumw2(kFALSE);
    
    // Fill individual channel histogram
    if (mMapHistAmpPerChannel.find(chID) != mMapHistAmpPerChannel.end() && 
        mMapHistAmpPerChannel[chID] != nullptr) {
      mMapHistAmpPerChannel[chID]->Add(projChannel.get());
    } else {
      ILOG(Warning, Support) << "Channel histogram for chID " << chID << " not found" << ENDM;
      continue;
    }
    
    // Add to combined histogram
    mHistAmpAll->Add(projChannel.get());
    
    // Add to normalized histogram with equal weight per channel
    const double integral = projChannel->Integral();
    if (integral > 0.0) {
      mHistAmpNormPerChannel->Add(projChannel.get(), 1.0 / integral);
    }
  }
  
  // Extract and set timestamp metadata
  long long timestamp = t.timestamp;
  
  if (mo) {
    for (const auto& metainfo : mo->getMetadataMap()) {
      if (metainfo.first == mTimestampMetaField) {
        try {
          timestamp = std::stoll(metainfo.second);
          ILOG(Debug, Support) << "Using timestamp from metadata: " << timestamp << ENDM;
        } catch (const std::exception& e) {
          ILOG(Warning, Support) << "Failed to parse timestamp from metadata: " << e.what() << ENDM;
        }
        break;
      }
    }
  }
  
  setTimestampToMOs(timestamp);
  
  ILOG(Info, Support) << "Successfully updated amplitude histograms for " << sNCHANNELS_PM << " channels" << ENDM;
}

void AmplitudePostProcTask::setTimestampToMOs(long long timestamp)
{
  try {
    const int nObjects = getObjectsManager()->getNumberPublishedObjects();
    for (int iObj = 0; iObj < nObjects; iObj++) {
      auto mo = getObjectsManager()->getMonitorObject(iObj);
      if (mo) {
        mo->addOrUpdateMetadata(mTimestampMetaField, std::to_string(timestamp));
      }
    }
    ILOG(Debug, Support) << "Set timestamp " << timestamp << " to " << nObjects << " monitor objects" << ENDM;
  } catch (const std::exception& e) {
    ILOG(Error, Support) << "Error setting timestamp to monitor objects: " << e.what() << ENDM;
  }
}

void AmplitudePostProcTask::finalize(Trigger, framework::ServiceRegistryRef)
{
  ILOG(Info, Support) << "Finalizing AmplitudePostProcTask" << ENDM;
}

} // namespace o2::quality_control_modules::fv0