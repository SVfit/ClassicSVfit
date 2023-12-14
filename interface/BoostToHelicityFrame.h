#ifndef TauAnalysis_ClassicSVfit_BoostToHelicityFrame_h
#define TauAnalysis_ClassicSVfit_BoostToHelicityFrame_h

#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"   // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // LorentzVector, Vector

#include <Math/Boost.h>                                           // ROOT::Math::Boost

class BoostToHelicityFrame
{
 public:
  BoostToHelicityFrame();
  ~BoostToHelicityFrame();

  void
  setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2);

  enum { kTauPlus, kTauMinus };
  classic_svFit::LorentzVector
  operator()(const classic_svFit::LorentzVector& p4, int tau);

 private:
  LorentzVector beamP4_;

  ROOT::Math::Boost boost_ttrf_;
  ROOT::Math::Boost boost_tprf_;
  ROOT::Math::Boost boost_tmrf_;

  classic_svFit::Vector r_;
  classic_svFit::Vector n_;
  classic_svFit::Vector k_;
};

#endif
