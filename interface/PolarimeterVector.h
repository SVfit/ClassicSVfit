#ifndef TauAnalysis_ClassicSVfit_PolarimeterVector_h
#define TauAnalysis_ClassicSVfit_PolarimeterVector_h

#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h"     // BoostToHelicityFrame
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"          // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"        // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoOneProng0Pi0.h"   // PolVecAlgoOneProng0Pi0
#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoOneProng1Pi0.h"   // PolVecAlgoOneProng1Pi0
#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoThreeProng0Pi0.h" // PolVecAlgoThreeProng0Pi0
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"        // Vector

namespace classic_svFit
{
  class PolarimeterVector
  {
   public:
    PolarimeterVector();
    ~PolarimeterVector();

    Vector
    operator()(const MeasuredTauLepton& measuredTauLepton, const FittedTauLepton& fittedTauLepton, int tau, 
               const BoostToHelicityFrame& boostToHelicityFrame) const;

   private:
    PolVecAlgoOneProng0Pi0 algoOneProng0Pi0_;
    PolVecAlgoOneProng1Pi0 algoOneProng1Pi0_;
    PolVecAlgoThreeProng0Pi0 algoThreeProng0Pi0_;
  };
}

#endif
