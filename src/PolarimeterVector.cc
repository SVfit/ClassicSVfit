#include "TauAnalysis/ClassicSVfit/interface/PolarimeterVector.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h" // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

namespace classic_svFit
{

PolarimeterVector::PolarimeterVector()
{}

PolarimeterVector::~PolarimeterVector()
{}

Vector
PolarimeterVector::operator()(const MeasuredTauLepton& measuredTauLepton, const FittedTauLepton& fittedTauLepton, int tau, 
                              const BoostToHelicityFrame& boostToHelicityFrame) const
{
  int tau_decayMode = measuredTauLepton.decayMode();
  Vector h;
  if ( tau_decayMode == reco::PFTau::kOneProng0PiZero )
  {
    h = algoOneProng0Pi0_(measuredTauLepton, tau, boostToHelicityFrame);
  }
  else if ( tau_decayMode == reco::PFTau::kOneProng1PiZero )
  {
    h = algoOneProng1Pi0_(measuredTauLepton, fittedTauLepton, tau, boostToHelicityFrame);
  }
  else if ( tau_decayMode == reco::PFTau::kThreeProng0PiZero )
  {
    h = algoThreeProng0Pi0_(measuredTauLepton, tau, boostToHelicityFrame);
  }
  return h;
}

}
