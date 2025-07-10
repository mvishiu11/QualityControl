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
/// \file   AmplitudeDistributionCheck.h
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Checks amplitude distributions for reasonable data content and statistics
///

#ifndef QC_MODULE_FV0_AMPLITUDEDISTRIBUTIONCHECK_H
#define QC_MODULE_FV0_AMPLITUDEDISTRIBUTIONCHECK_H

#include "QualityControl/CheckInterface.h"
#include "FV0Base/Constants.h"
#include <TF1.h>

namespace o2::quality_control_modules::fv0
{

/// \brief Checks amplitude distribution histograms for data quality and completeness
/// 
/// This check analyzes amplitude histograms to ensure:
/// 1. Sufficient statistics for meaningful analysis
/// 2. Reasonable distribution shapes (not all zeros, not dominated by overflow/underflow)
/// 3. Consistent behavior across different amplitude summaries
/// 4. Detection of potential hardware or readout issues
///
/// This check consumes multiple histogram types: AmpAllChannels, AmpNormPerChannel,
/// and individual channel histograms, ensuring pipeline completion while providing
/// comprehensive quality assessment.
///
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
class AmplitudeDistributionCheck : public o2::quality_control::checker::CheckInterface
{
public:
    AmplitudeDistributionCheck() = default;
    ~AmplitudeDistributionCheck() override = default;

    void configure() override;
    Quality check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap) override;
    void beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult = Quality::Null) override;
    std::string getAcceptedType() override;

    ClassDefOverride(AmplitudeDistributionCheck, 1);

private:
    // Configuration parameters
    int mMinEntries{1000};                   ///< Minimum entries required for good quality
    int mMinEntriesWarning{100};             ///< Minimum entries for warning quality
    float mMaxEmptyChannelFraction{0.1};     ///< Maximum fraction of empty channels allowed
    float mMaxOverflowFraction{0.05};        ///< Maximum fraction in overflow bin
    float mMaxUnderflowFraction{0.05};       ///< Maximum fraction in underflow bin
    float mMinMeanAmplitude{10.0};           ///< Minimum expected mean amplitude
    float mMaxMeanAmplitude{1000.0};         ///< Maximum expected mean amplitude
    bool mCheckIndividualChannels{true};     ///< Whether to check individual channel histograms
    bool mEnableMessage{true};               ///< Whether to add quality messages to plots
    int mMaxChannelsToAnalyze{48};           ///< Maximum number of individual channels to analyze

    // Analysis results
    struct AnalysisResults {
        int totalHistogramsAnalyzed{0};
        int histogramsWithSufficientStats{0};
        int emptyChannels{0};
        int channelsWithOverflow{0};
        int channelsWithUnderflow{0};
        double totalEntries{0};
        double meanAmplitudeOverall{0.0};
        double minChannelMean{999999.0};
        double maxChannelMean{0.0};
        bool foundAmpAllChannels{false};
        bool foundAmpNormPerChannel{false};
        std::string firstProblemChannel{""};
    } mResults;

    // Helper methods for different types of analysis
    Quality analyzeMainHistograms(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap);
    Quality analyzeIndividualChannels(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap);
    Quality combineQualities(Quality main, Quality individual) const;
    bool isHistogramReasonable(TH1* hist, const std::string& name) const;
    void updateResultsFromHistogram(TH1* hist, const std::string& name);
};

} // namespace o2::quality_control_modules::fv0

#endif // QC_MODULE_FV0_AMPLITUDEDISTRIBUTIONCHECK_H