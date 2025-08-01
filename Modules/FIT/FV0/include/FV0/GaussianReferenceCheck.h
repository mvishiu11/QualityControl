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
/// \file   GaussianReferenceCheck.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Checks if Gaussian mean values are consistent with reference line at 1.0
///

#ifndef QC_MODULE_FV0_GAUSSIANREFERENCECHECK_H
#define QC_MODULE_FV0_GAUSSIANREFERENCECHECK_H

#include "QualityControl/CheckInterface.h"
#include "FV0Base/Constants.h"

namespace o2::quality_control_modules::fv0
{

/// \brief Checks if scaled Gaussian means are within acceptable range of reference value
/// 
/// This check analyzes the GaussianSummary/MeanVsChannel graph to ensure that:
/// 1. Error bars of measurement points overlap with the reference line at y=1.0
/// 2. The majority of channels show consistent response
/// 3. No systematic deviations that might indicate calibration or hardware issues
///
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
class GaussianReferenceCheck : public o2::quality_control::checker::CheckInterface
{
public:
    GaussianReferenceCheck() = default;
    ~GaussianReferenceCheck() override = default;

    void configure() override;
    Quality check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap) override;
    void beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult = Quality::Null) override;
    std::string getAcceptedType() override;

    ClassDefOverride(GaussianReferenceCheck, 1);

private:
    // Configuration parameters for quality thresholds
    float mReferenceValue{1.0};           ///< Expected reference value for scaled means
    float mToleranceWarning{0.2};         ///< Warning threshold for deviation from reference
    float mToleranceError{0.4};           ///< Error threshold for deviation from reference  
    float mMinChannelFractionGood{0.8};   ///< Minimum fraction of channels that must be good
    float mMinChannelFractionMedium{0.6}; ///< Minimum fraction of channels for medium quality
    bool mRequireErrorBarOverlap{true};   ///< Whether to require error bar overlap with reference
    bool mEnableMessage{true};            ///< Whether to add quality message to plots
    
    // Analysis results
    int mTotalChannelsAnalyzed{0};
    int mChannelsWithinReference{0};
    int mChannelsWithErrorBarOverlap{0};
    float mMeanDeviation{0.0};
    float mMaxDeviation{0.0};
    
    // Helper methods for analysis
    bool isPointWithinTolerance(double value, double error, double tolerance) const;
    bool doesErrorBarOverlapReference(double value, double error) const;
};

} // namespace o2::quality_control_modules::fv0

#endif // QC_MODULE_FV0_GAUSSIANREFERENCECHECK_H