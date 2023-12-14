#ifndef TauAnalysis_ClassicSVfit_PolVecAlgoOneProng0Pi0_h
#define TauAnalysis_ClassicSVfit_PolVecAlgoOneProng0Pi0_h

#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h" // BoostToHelicityFrame
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLeptonLT.h"  // MeasuredTauLeptonLT
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"    // Vector

class PolVecAlgoOneProng0Pi0
{
 public:
  PolVecAlgoOneProng0Pi0();
  ~PolVecAlgoOneProng0Pi0();

  classic_svFit::Vector
  operator()(const MeasuredTauLeptonLT& measuredTauLepton, int tau,
             const BoostToHelicityFrame& boostToHelicityFrame) const;
};

#endif
