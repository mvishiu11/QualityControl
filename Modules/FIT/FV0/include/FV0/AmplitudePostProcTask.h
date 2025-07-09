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
#include "QualityControl/DatabaseInterface.h"
#include "FV0Base/Constants.h"
#include "CCDB/CcdbApi.h"

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

    class AmplitudePostProcTask final : public quality_control::postprocessing::PostProcessingInterface
    {
        public:
            AmplitudePostProcTask()  = default;
            ~AmplitudePostProcTask() override = default;

            // QC framework hooks
            void configure  (const boost::property_tree::ptree& config) override;
            void initialize (quality_control::postprocessing::Trigger trigger,
                            framework::ServiceRegistryRef services)   override;
            void update     (quality_control::postprocessing::Trigger trigger,
                            framework::ServiceRegistryRef services)   override;
            void finalize   (quality_control::postprocessing::Trigger trigger,
                            framework::ServiceRegistryRef services)   override;

        private:
            // helper utilities
            void reset();
            void setTimestampToMOs(long long timestamp);

            //  Configuration & constants
            static constexpr std::size_t sNCHANNELS_PM =
                o2::fv0::Constants::nFv0ChannelsPlusRef;

            std::string mPathDigitQcTask;                   ///< Source MO path
            std::string mCcdbUrl;                           ///< CCDB URL
            std::string mTimestampMetaField{"timestampTF"}; ///< Metadata key
            int    mAmpMin  {-100};                         ///< Histogram ADC min
            int    mAmpMax  {4100};                         ///< Histogram ADC max
            int    mAmpBins {4200};                         ///< Histogram bins
            double mSliceFrac{0.25};                        ///< Fit window half-width (fraction of peak)

            //  QC managed objects
            o2::quality_control::repository::DatabaseInterface* mDatabase{nullptr};
            o2::ccdb::CcdbApi mCcdbApi;

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
    };

}   // namespace o2::quality_control_modules::fv0
#endif // QC_MODULE_FV0_AMPLITUDEPOSTPROCTASK_H