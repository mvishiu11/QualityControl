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
/// \file   AgingLaserPostProcTask.h
/// \author Andreas Molander <andreas.molander@cern.ch>, Jakub Muszyński <jakub.milosz.muszynski@cern.ch>
///

#ifndef QC_MODULE_FT0_AGINGLASERPOSTPROC_H
#define QC_MODULE_FT0_AGINGLASERPOSTPROC_H

#include "QualityControl/PostProcessingInterface.h"
#include <FT0Base/Constants.h>

#include <TH1F.h>
#include <memory>
#include <vector>

namespace o2::quality_control_modules::ft0
{

/// \brief Post-processing task that derives a per-channel
///        **weighted-mean amplitude, normalised to reference channels**.
///        The algorithm:
///        1. For every configured reference channel:
///             * project the corresponding slice of `AmpPerChannel` (TH2)
///             * in ADC range `[adcSearchMin, adcSearchMax]` find the max bin centre *x_max*
///             * fit a Gaussian in `[ (1-a)*x_max , (1+b)*x_max ]` → μ<sub>ref</sub>
///        2. `norm = average( μ<sub>ref</sub> )`
///        3. For **every** channel (detector + reference):
///             * find global *x_max* of the slice
///             * compute a **weighted mean** inside the same fractional window
///             * value = (weighted mean) / norm
///        4. Store all values in one TH1F (bin = channel).
class AgingLaserPostProcTask final : public quality_control::postprocessing::PostProcessingInterface
{
 public:
  AgingLaserPostProcTask() = default;
  ~AgingLaserPostProcTask() override;

  void initialize(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;
  void update(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;
  void finalize(quality_control::postprocessing::Trigger, framework::ServiceRegistryRef) override;

 private:
  /// Maximum number of FT0 photo-multiplier channels
  static constexpr std::size_t sNCHANNELS_PM = o2::ft0::Constants::sNCHANNELS_PM;

  /* ---------------- task parameters (configurable) ---------------- */
  std::vector<uint8_t> mDetectorChIDs;  ///< Detector-side channels
  std::vector<uint8_t> mReferenceChIDs; ///< Reference (laser monitor) channels

  double mADCSearchMin = 150.; ///< lower edge of peak-search window (ADC)
  double mADCSearchMax = 600.; ///< upper edge of peak-search window (ADC)
  double mFracWindowA = 0.25;  ///< low fractional window parameter *a*
  double mFracWindowB = 0.25;  ///< high fractional window parameter *b*

  /* ---------------- output histogram ---------------- */
  std::unique_ptr<TH1F> mAmpVsChNormWeightedMeanA;
  std::unique_ptr<TH1F> mAmpVsChNormWeightedMeanC;
};

} // namespace o2::quality_control_modules::ft0

#endif // QC_MODULE_FT0_AGINGLASERPOSTPROC_H
