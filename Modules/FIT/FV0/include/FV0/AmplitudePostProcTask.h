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

#include "QualityControl/PostProcessingInterface.h"
#include "QualityControl/DatabaseInterface.h"
#include "FV0Base/Constants.h"
#include "CCDB/CcdbApi.h"

#include <memory>
#include <string>
#include <map>

class TH1F;
class TH2F;

namespace o2::quality_control_modules::fv0
{

/// \brief Post-processing task for FV0 amplitude analysis
/// \details Creates individual amplitude histograms for each channel and combined histograms
///          for overall analysis
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
class AmplitudePostProcTask final : public quality_control::postprocessing::PostProcessingInterface
{
 public:
  /// \brief Constructor
  AmplitudePostProcTask() = default;
  
  /// \brief Destructor
  ~AmplitudePostProcTask() override;

  /// \brief Configuration of the task
  void configure(const boost::property_tree::ptree& config) override;
  
  /// \brief Initialization of the task
  void initialize(quality_control::postprocessing::Trigger trigger, framework::ServiceRegistryRef services) override;
  
  /// \brief Update of the task, triggered by external events
  void update(quality_control::postprocessing::Trigger trigger, framework::ServiceRegistryRef services) override;
  
  /// \brief Finalization of the task
  void finalize(quality_control::postprocessing::Trigger trigger, framework::ServiceRegistryRef services) override;

 private:
  /// \brief Reset all histogram objects
  void reset();
  
  /// \brief Set timestamp metadata to all monitor objects
  /// \param timestamp Unix timestamp in milliseconds
  void setTimestampToMOs(long long timestamp);

  // Constants
  static constexpr std::size_t sNCHANNELS_PM = o2::fv0::Constants::nFv0ChannelsPlusRef;
  
  // Configuration parameters
  std::string mPathDigitQcTask;                      ///< Path to the digit QC task
  std::string mCcdbUrl;                              ///< CCDB URL
  std::string mTimestampMetaField{"timestampTF"};    ///< Metadata field for timestamp
  int mAmpMin{-100};                                 ///< Minimum amplitude for histograms
  int mAmpMax{4100};                                 ///< Maximum amplitude for histograms  
  int mAmpBins{4200};                                ///< Number of bins for amplitude histograms

  // Database interfaces
  o2::quality_control::repository::DatabaseInterface* mDatabase = nullptr; ///< Database interface
  o2::ccdb::CcdbApi mCcdbApi;                                               ///< CCDB API

  // Monitor objects
  std::map<unsigned int, std::unique_ptr<TH1F>> mMapHistAmpPerChannel; ///< Individual channel amplitude histograms
  std::unique_ptr<TH1F> mHistAmpAll;                                   ///< Combined amplitude histogram for all channels
  std::unique_ptr<TH1F> mHistAmpNormPerChannel;                        ///< Normalized amplitude histogram per channel
};

} // namespace o2::quality_control_modules::fv0

#endif // QC_MODULE_FV0_AMPLITUDEPOSTPROCTASK_H