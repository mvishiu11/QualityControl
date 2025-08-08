#ifndef QC_MODULE_FV0_FV0MIPCHECK_H
#define QC_MODULE_FV0_FV0MIPCHECK_H

#include "QualityControl/CheckInterface.h"
#include "CommonConstants/LHCConstants.h"

namespace o2::quality_control_modules::fv0
{

class FV0MIPCheck : public o2::quality_control::checker::CheckInterface
{
 public:
  FV0MIPCheck()  = default;
  ~FV0MIPCheck() override = default;

  /* CheckInterface --------------------------------------------------- */
  void configure() override;
  Quality check(std::map<std::string,
        std::shared_ptr<MonitorObject>>* moMap) override;
  void beautify(std::shared_ptr<MonitorObject> mo,
                Quality res = Quality::Null) override;
  std::string getAcceptedType() override { return "TGraphErrors"; }

 private:
  double mReferenceADC {15.};   ///< expected gain
  double mToleranceADC {1.};    ///< |y-ref| that triggers warning
  bool   mDrawWindow   {true};  ///< paint the band?

  double mMaxDeviation {-1.};   ///< filled in check(), shown in label
};

} // namespace o2::quality_control_modules::fv0

#endif
