//  FV0MIPCheck.cxx
#include "FV0/FV0MIPCheck.h"

#include "QualityControl/MonitorObject.h"
#include "QualityControl/Quality.h"
#include "QualityControl/QcInfoLogger.h"
#include "Common/Utils.h"

#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TAxis.h> 
#include <algorithm>

using namespace o2::quality_control;

namespace o2::quality_control_modules::fv0
{

void FV0MIPCheck::configure()
{
  mReferenceADC = std::stod(mCustomParameters.atOrDefaultValue("referenceADC", "15"));
  mToleranceADC = std::stod(mCustomParameters.atOrDefaultValue("toleranceADC", "1"));
  mDrawWindow   = o2::quality_control_modules::common::getFromConfig(
                    mCustomParameters, "drawWindow", true);
}

Quality FV0MIPCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  Quality res = Quality::Null;
  TGraphErrors* graph = nullptr;
  for (auto& [name, mo] : *moMap) {
    graph = dynamic_cast<TGraphErrors*>(mo->getObject());
    if (graph) break;
  }
  if (!graph) {
    ILOG(Warning, Support) << "FV0MIPCheck: no TGraphErrors in moMap - returning Null" << ENDM;
    return res;
  }

  res = Quality::Good;
  const int n = graph->GetN();
  mMaxDeviation = 0.;

  for (int i = 0; i < n; ++i) {
    const double y   = graph->GetPointY(i);
    const double dev = std::abs(y - mReferenceADC);
    mMaxDeviation = std::max(mMaxDeviation, dev);
    if (dev > mToleranceADC) { res.set(Quality::Medium); }
  }

  res.addMetadata("maxDeviationADC", std::to_string(mMaxDeviation));
  return res;
}

void FV0MIPCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality q)
{
  auto* g = dynamic_cast<TGraphErrors*>(mo->getObject());
  if (!g) return;

  if (mDrawWindow) {
    const double xmin = g->GetXaxis()->GetXmin();
    const double xmax = g->GetXaxis()->GetXmax();

    auto* up = new TLine(xmin,
                         mReferenceADC + mToleranceADC,
                         xmax,
                         mReferenceADC + mToleranceADC);

    auto* dn = new TLine(xmin,
                         mReferenceADC - mToleranceADC,
                         xmax,
                         mReferenceADC - mToleranceADC);

    up->SetLineColor(kOrange);
    dn->SetLineColor(kOrange);
    up->SetLineStyle(kDashed);
    dn->SetLineStyle(kDashed);

    g->GetListOfFunctions()->Add(up);
    g->GetListOfFunctions()->Add(dn);
    g->GetYaxis()->SetRangeUser(mReferenceADC - 2, mReferenceADC + 2);
  }

  TPaveText* msg = new TPaveText(0.15, 0.82, 0.85, 0.92, "NDC");
  msg->SetTextAlign(12);
  msg->AddText(Form("deviation = %.2f  (tolerance = %.2f)\n", mMaxDeviation, mToleranceADC));
  if (q != Quality::Good) {
    msg->AddText("FV0 gain calibration failing, please make log entry and tag FV0, FIT");
  } else {
    msg->AddText("FV0 gain calibration OK");
  }
  msg->SetFillColor(q == Quality::Good   ? kGreen :
                    q == Quality::Medium ? kOrange : kGray);
  g->GetListOfFunctions()->Add(msg);
}

} // namespace o2::quality_control_modules::fv0
