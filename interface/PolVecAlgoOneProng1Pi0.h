#ifndef TauAnalysis_ClassicSVfit_PolVecAlgoOneProng1Pi0_h
#define TauAnalysis_ClassicSVfit_PolVecAlgoOneProng1Pi0_h

#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h" // BoostToHelicityFrame
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"    // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"      // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"    // Vector

namespace classic_svFit
{
  class PolVecAlgoOneProng1Pi0
  {
   public:
    PolVecAlgoOneProng1Pi0();
    ~PolVecAlgoOneProng1Pi0();

    Vector
    operator()(const MeasuredTauLepton& measuredTauLepton, const FittedTauLepton& fittedTauLepton, int tau, 
               const BoostToHelicityFrame& boostToHelicityFrame) const;
  };
}

#endif
