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
/// \file   AmplitudeDistributionCheck.cxx
/// \author Jakub Muszyński jakub.milosz.muszynski@cern.ch
/// \brief  Checks amplitude distributions for reasonable data content and statistics
///

#include "FV0/AmplitudeDistributionCheck.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/Quality.h"
#include "QualityControl/QcInfoLogger.h"
#include "Common/Utils.h"

// ROOT includes
#include <TH1.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TList.h>
#include <TMath.h>

#include <DataFormatsQualityControl/FlagType.h>
#include <DataFormatsQualityControl/FlagTypeFactory.h>

using namespace std;
using namespace o2::quality_control;

namespace o2::quality_control_modules::fv0
{

void AmplitudeDistributionCheck::configure()
{
    // Configure minimum entries thresholds
    if (auto param = mCustomParameters.find("minEntries"); param != mCustomParameters.end()) {
        mMinEntries = stoi(param->second);
        ILOG(Debug, Support) << "configure() : using minEntries = " << mMinEntries << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default minEntries = " << mMinEntries << ENDM;
    }

    if (auto param = mCustomParameters.find("minEntriesWarning"); param != mCustomParameters.end()) {
        mMinEntriesWarning = stoi(param->second);
        ILOG(Debug, Support) << "configure() : using minEntriesWarning = " << mMinEntriesWarning << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default minEntriesWarning = " << mMinEntriesWarning << ENDM;
    }

    // Configure quality thresholds for distribution analysis
    if (auto param = mCustomParameters.find("maxEmptyChannelFraction"); param != mCustomParameters.end()) {
        mMaxEmptyChannelFraction = stof(param->second);
        ILOG(Debug, Support) << "configure() : using maxEmptyChannelFraction = " << mMaxEmptyChannelFraction << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default maxEmptyChannelFraction = " << mMaxEmptyChannelFraction << ENDM;
    }

    if (auto param = mCustomParameters.find("maxOverflowFraction"); param != mCustomParameters.end()) {
        mMaxOverflowFraction = stof(param->second);
        ILOG(Debug, Support) << "configure() : using maxOverflowFraction = " << mMaxOverflowFraction << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default maxOverflowFraction = " << mMaxOverflowFraction << ENDM;
    }

    // Configure amplitude range expectations
    if (auto param = mCustomParameters.find("minMeanAmplitude"); param != mCustomParameters.end()) {
        mMinMeanAmplitude = stof(param->second);
        ILOG(Debug, Support) << "configure() : using minMeanAmplitude = " << mMinMeanAmplitude << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default minMeanAmplitude = " << mMinMeanAmplitude << ENDM;
    }

    if (auto param = mCustomParameters.find("maxMeanAmplitude"); param != mCustomParameters.end()) {
        mMaxMeanAmplitude = stof(param->second);
        ILOG(Debug, Support) << "configure() : using maxMeanAmplitude = " << mMaxMeanAmplitude << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default maxMeanAmplitude = " << mMaxMeanAmplitude << ENDM;
    }

    // Configure analysis scope and messaging
    mCheckIndividualChannels = o2::quality_control_modules::common::getFromConfig(mCustomParameters, "checkIndividualChannels", true);
    mEnableMessage = o2::quality_control_modules::common::getFromConfig(mCustomParameters, "enableMessage", true);

    if (auto param = mCustomParameters.find("maxChannelsToAnalyze"); param != mCustomParameters.end()) {
        mMaxChannelsToAnalyze = stoi(param->second);
        ILOG(Debug, Support) << "configure() : using maxChannelsToAnalyze = " << mMaxChannelsToAnalyze << ENDM;
    } else {
        ILOG(Debug, Support) << "configure() : using default maxChannelsToAnalyze = " << mMaxChannelsToAnalyze << ENDM;
    }
}

Quality AmplitudeDistributionCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
    // Reset analysis results for this check cycle
    mResults = AnalysisResults();
    
    ILOG(Debug, Support) << "Starting amplitude distribution check with " << moMap->size() << " monitoring objects" << ENDM;

    // Analyze the main summary histograms first
    Quality mainQuality = analyzeMainHistograms(moMap);
    
    // Analyze individual channel histograms if enabled and available
    Quality individualQuality = Quality::Good; 
    if (mCheckIndividualChannels) {
        individualQuality = analyzeIndividualChannels(moMap);
    }

    // Combine the two quality assessments
    Quality finalQuality = combineQualities(mainQuality, individualQuality);

    // Add metadata for trending and debugging
    finalQuality.addMetadata("totalHistogramsAnalyzed", std::to_string(mResults.totalHistogramsAnalyzed));
    finalQuality.addMetadata("histogramsWithSufficientStats", std::to_string(mResults.histogramsWithSufficientStats));
    finalQuality.addMetadata("emptyChannels", std::to_string(mResults.emptyChannels));
    finalQuality.addMetadata("totalEntries", std::to_string(static_cast<int>(mResults.totalEntries)));
    finalQuality.addMetadata("meanAmplitudeOverall", std::to_string(mResults.meanAmplitudeOverall));
    finalQuality.addMetadata("foundAmpAllChannels", mResults.foundAmpAllChannels ? "true" : "false");
    finalQuality.addMetadata("foundAmpNormPerChannel", mResults.foundAmpNormPerChannel ? "true" : "false");

    ILOG(Info, Support) << "Amplitude distribution check completed: analyzed " 
                       << mResults.totalHistogramsAnalyzed << " histograms, "
                       << mResults.histogramsWithSufficientStats << " with sufficient statistics, "
                       << "quality=" << finalQuality << ENDM;

    return finalQuality;
}

Quality AmplitudeDistributionCheck::analyzeMainHistograms(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
    Quality result = Quality::Good;
    std::vector<std::string> problemDescriptions;

    // Main amplitude histograms that the post-processing task creates
    for (auto& [moName, mo] : *moMap) {
        TH1* hist = dynamic_cast<TH1*>(mo->getObject());
        if (!hist) continue;

        std::string name = mo->getName();
        
        // Check if this is one of main amplitude histograms
        if (name.find("AmpAllChannels") != std::string::npos) {
            mResults.foundAmpAllChannels = true;
            mResults.totalHistogramsAnalyzed++;
            
            if (!isHistogramReasonable(hist, name)) {
                problemDescriptions.push_back("AmpAllChannels has insufficient or problematic data");
                result = Quality::Bad;
            } else {
                mResults.histogramsWithSufficientStats++;
            }
            updateResultsFromHistogram(hist, name);
            
        } else if (name.find("AmpNormPerChannel") != std::string::npos) {
            mResults.foundAmpNormPerChannel = true;
            mResults.totalHistogramsAnalyzed++;
            
            if (!isHistogramReasonable(hist, name)) {
                problemDescriptions.push_back("AmpNormPerChannel has insufficient or problematic data");
                if (result == Quality::Good) result = Quality::Medium; // Less critical than AmpAllChannels
            } else {
                mResults.histogramsWithSufficientStats++;
            }
            updateResultsFromHistogram(hist, name);
        }
    }

    // Check if we found the expected main histograms
    if (!mResults.foundAmpAllChannels) {
        problemDescriptions.push_back("Missing expected AmpAllChannels histogram");
        result = Quality::Bad;
    }
    
    if (!mResults.foundAmpNormPerChannel) {
        problemDescriptions.push_back("Missing expected AmpNormPerChannel histogram");
        if (result == Quality::Good) result = Quality::Medium;
    }

    // Add flags for any problems found
    for (const auto& problem : problemDescriptions) {
        result.addFlag(FlagTypeFactory::Unknown(), problem);
    }

    return result;
}

Quality AmplitudeDistributionCheck::analyzeIndividualChannels(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
    Quality result = Quality::Good;
    int channelsAnalyzed = 0;
    int problematicChannels = 0;

    // Individual channel amplitude histograms
    for (auto& [moName, mo] : *moMap) {
        // Limit the number of individual channels we analyze to avoid overwhelming the check
        if (channelsAnalyzed >= mMaxChannelsToAnalyze) break;

        TH1* hist = dynamic_cast<TH1*>(mo->getObject());
        if (!hist) continue;

        std::string name = mo->getName();
        
        // Check if this is an individual channel histogram
        if (name.find("AmplitudePerChannel") != std::string::npos && name.find("Amp_ch") != std::string::npos) {
            channelsAnalyzed++;
            mResults.totalHistogramsAnalyzed++;
            
            if (hist->GetEntries() == 0) {
                mResults.emptyChannels++;
            } else if (!isHistogramReasonable(hist, name)) {
                problematicChannels++;
                if (mResults.firstProblemChannel.empty()) {
                    mResults.firstProblemChannel = name;
                }
            } else {
                mResults.histogramsWithSufficientStats++;
            }
            
            updateResultsFromHistogram(hist, name);
        }
    }

    // Assess quality based on the fraction of problematic channels
    if (channelsAnalyzed > 0) {
        float emptyChannelFraction = static_cast<float>(mResults.emptyChannels) / channelsAnalyzed;
        float problematicChannelFraction = static_cast<float>(problematicChannels) / channelsAnalyzed;

        if (emptyChannelFraction > mMaxEmptyChannelFraction) {
            result = Quality::Bad;
            result.addFlag(FlagTypeFactory::Unknown(), 
                          Form("Too many empty channels: %.1f%% (threshold: %.1f%%)", 
                               emptyChannelFraction * 100, mMaxEmptyChannelFraction * 100));
        } else if (problematicChannelFraction > 0.2) {
            result = Quality::Medium;
            result.addFlag(FlagTypeFactory::Unknown(),
                          Form("%.1f%% of individual channels show data quality issues", 
                               problematicChannelFraction * 100));
        }
    }

    ILOG(Debug, Support) << "Individual channel analysis: " << channelsAnalyzed << " channels analyzed, "
                        << mResults.emptyChannels << " empty, " << problematicChannels << " problematic" << ENDM;

    return result;
}

std::string AmplitudeDistributionCheck::getAcceptedType() 
{ 
    return "TH1"; 
}

bool AmplitudeDistributionCheck::isHistogramReasonable(TH1* hist, const std::string& name) const
{
    if (!hist) return false;

    // Check for sufficient statistics
    if (hist->GetEntries() < mMinEntriesWarning) {
        ILOG(Debug, Support) << "Histogram " << name << " has insufficient entries: " << hist->GetEntries() << ENDM;
        return false;
    }

    // Check for reasonable mean
    double mean = hist->GetMean();
    if (mean < mMinMeanAmplitude || mean > mMaxMeanAmplitude) {
        ILOG(Debug, Support) << "Histogram " << name << " has unreasonable mean: " << mean << ENDM;
        return false;
    }

    // Check overflow and underflow fractions
    double totalEntries = hist->GetEntries();
    if (totalEntries > 0) {
        double overflowFraction = hist->GetBinContent(hist->GetNbinsX() + 1) / totalEntries;
        double underflowFraction = hist->GetBinContent(0) / totalEntries;
        
        if (overflowFraction > mMaxOverflowFraction || underflowFraction > mMaxUnderflowFraction) {
            ILOG(Debug, Support) << "Histogram " << name << " has excessive overflow/underflow: " 
                                << overflowFraction << "/" << underflowFraction << ENDM;
            return false;
        }
    }

    return true;
}

void AmplitudeDistributionCheck::updateResultsFromHistogram(TH1* hist, const std::string& name)
{
    if (!hist) return;

    mResults.totalEntries += hist->GetEntries();
    double mean = hist->GetMean();
    
    // Update overall statistics
    if (mean < mResults.minChannelMean) mResults.minChannelMean = mean;
    if (mean > mResults.maxChannelMean) mResults.maxChannelMean = mean;
    
    // Calculate running average of mean amplitudes
    static double totalMeanSum = 0.0;
    static int meanCount = 0;
    totalMeanSum += mean;
    meanCount++;
    mResults.meanAmplitudeOverall = totalMeanSum / meanCount;

    // Check for overflow/underflow issues
    if (hist->GetEntries() > 0) {
        double overflowFraction = hist->GetBinContent(hist->GetNbinsX() + 1) / hist->GetEntries();
        double underflowFraction = hist->GetBinContent(0) / hist->GetEntries();
        
        if (overflowFraction > mMaxOverflowFraction) mResults.channelsWithOverflow++;
        if (underflowFraction > mMaxUnderflowFraction) mResults.channelsWithUnderflow++;
    }
}

Quality AmplitudeDistributionCheck::combineQualities(Quality main, Quality individual) const
{
    // The main histograms are more critical than individual channels
    // Return the worse of the two, with main quality taking precedence in ties
    if (main == Quality::Bad || individual == Quality::Bad) {
        return Quality::Bad;
    } else if (main == Quality::Medium || individual == Quality::Medium) {
        return Quality::Medium;
    } else if (main == Quality::Null || individual == Quality::Null) {
        return Quality::Null;
    } else {
        return Quality::Good;
    }
}

void AmplitudeDistributionCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
    auto* hist = dynamic_cast<TH1*>(mo->getObject());
    if (hist == nullptr) {
        ILOG(Warning, Devel) << "Could not cast " << mo->getName() << " to TH1*, skipping beautification" << ENDM;
        return;
    }

    if (!mEnableMessage) {
        return;
    }

    // Only add messages to the main amplitude histograms to avoid clutter
    std::string name = mo->getName();
    if (name.find("AmpAllChannels") == std::string::npos && 
        name.find("AmpNormPerChannel") == std::string::npos) {
        return;
    }

    TPaveText* msg = new TPaveText(0.15, 0.75, 0.85, 0.95, "NDC");
    hist->GetListOfFunctions()->Add(msg);
    msg->SetName(Form("%s_msg", mo->GetName()));
    msg->Clear();
    msg->SetTextAlign(12);
    msg->SetTextSize(0.025);

    // Create informative summary message
    std::string summary = Form("Analyzed %d histograms, %d with sufficient stats, %d empty channels",
                              mResults.totalHistogramsAnalyzed, mResults.histogramsWithSufficientStats, 
                              mResults.emptyChannels);
    
    std::string details = Form("Total entries: %.0f, Mean amplitude: %.1f", 
                              mResults.totalEntries, mResults.meanAmplitudeOverall);

    // Set color and additional messages based on quality
    if (checkResult == Quality::Good) {
        msg->SetFillColor(kGreen);
        msg->AddText((summary + " >> Quality::Good <<").c_str());
        msg->AddText(details.c_str());
        msg->AddText("All amplitude distributions look reasonable");
    } else if (checkResult == Quality::Bad) {
        msg->SetFillColor(kRed);
        msg->AddText((summary + " >> Quality::Bad <<").c_str());
        msg->AddText(details.c_str());
        if (!mResults.firstProblemChannel.empty()) {
            msg->AddText(Form("First problematic channel: %s", mResults.firstProblemChannel.c_str()));
        }
        msg->AddText("Check detector readout and data acquisition!");
    } else if (checkResult == Quality::Medium) {
        msg->SetFillColor(kOrange);
        msg->AddText((summary + " >> Quality::Medium <<").c_str());
        msg->AddText(details.c_str());
        msg->AddText("Some amplitude distributions show minor issues");
    } else if (checkResult == Quality::Null) {
        msg->SetFillColor(kGray);
        msg->AddText("No valid amplitude data available >> Quality::Null <<");
    }
}

} // namespace o2::quality_control_modules::fv0