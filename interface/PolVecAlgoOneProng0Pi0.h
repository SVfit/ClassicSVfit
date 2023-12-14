#ifndef TauAnalysis_ClassicSVfit_PolVecAlgoOneProng0Pi0_h
#define TauAnalysis_ClassicSVfit_PolVecAlgoOneProng0Pi0_h

#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h" // BoostToHelicityFrame
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"    // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"    // Vector

namespace classic_svFit
{
  class PolVecAlgoOneProng0Pi0
  {
   public:
    PolVecAlgoOneProng0Pi0();
    ~PolVecAlgoOneProng0Pi0();

    Vector
    operator()(const MeasuredTauLepton& measuredTauLepton, int tau,
               const BoostToHelicityFrame& boostToHelicityFrame) const;
  };
}

#endif
