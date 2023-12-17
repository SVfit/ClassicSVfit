#ifndef TauAnalysis_ClassicSVfit_BoostToHelicityFrame_h
#define TauAnalysis_ClassicSVfit_BoostToHelicityFrame_h

#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"   // FittedTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // LorentzVector, Vector

#include <Math/Boost.h>                                           // ROOT::Math::Boost

namespace classic_svFit
{
  class BoostToHelicityFrame
  {
   public:
    BoostToHelicityFrame();
    ~BoostToHelicityFrame();

    void
    setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2);

    enum { kTauPlus, kTauMinus };
    LorentzVector
    operator()(const LorentzVector& p4, int tau) const;

    

   private:
    LorentzVector beamP4_;

    ROOT::Math::Boost boost_ttrf_;
    ROOT::Math::Boost boost_tprf_;
    ROOT::Math::Boost boost_tmrf_;

    Vector r_;
    Vector n_;
    Vector k_;
  };
}

#endif
