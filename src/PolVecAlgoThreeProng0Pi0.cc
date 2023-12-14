#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoOneProng1Pi0.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // fixNeutrinoMass(), gamma_va, tauLeptonMass

#include <assert.h>

PolVecAlgoThreeProng0Pi0::PolVecAlgoThreeProng0Pi0()
{}

PolVecAlgoThreeProng0Pi0::~PolVecAlgoThreeProng0Pi0()
{}

namespace
{
  bool
  isHigherPt(const MeasuredHadTauDecayProduct* particle1, const MeasuredHadTauDecayProduct* particle2)
  {
    return particle1->pt() > particle2->pt();
  }

  TLorentzVector
  convert_to_TLorentzVector(const LorentzVector& p4)
  {
    return TLorentzVector(p4.px(), p4.py(), p4.pz(), p4.energy());
  }
}

classic_svFit::Vector
PolVecAlgoOneProng1Pi0::operator()(const MeasuredTauLeptonLT& measuredTauLepton, const FittedTauLepton& fittedTauLepton, 
                                   const BoostToHelicityFrame& boostToHelicityFrame) const
{
  const std::vector<MeasuredHadTauDecayProduct>& daughters = measuredTauLepton.measuredHadTauDecayProducts();
  std::vector<const MeasuredHadTauDecayProduct*> chsPlus;
  std::vector<const MeasuredHadTauDecayProduct*> chsMinus;
  int charge_sum = 0;
  for ( const KinematicParticle& daughter : daughters )
  {
    if      ( daughter->charge() > 0 ) chsPlus.push_back(&daughter);
    else if ( daughter->charge() < 0 ) chsMinus.push_back(&daughter);
    charge_sum += daughter.charge();
  }
  if ( !((chsPlus.size() + chsMinus.size()) == 3 && std::abs(charge_sum) == 1) )
  {
    std::cerr << "ERROR: Failed to find three charged pions with |sum(charge)| = 1 !!" << std::endl;
    assert(0);
  }
  if ( charge_sum > 0 )
  {
    chsOS = chsMinus;
    chsSS = chsPlus;
  }
  else
  {
    chsOS = chsPlus;
    chsSS = chsMinus;
  }
  assert(chsOS.size() == 1 && chsSS.size() == 2);
  const LorentzVector& chOSP4 = chsOS[0]->p4();
  std::sort(chsSS.begin(), chsSS.end(), isHigherPt);
  const LorentzVector& chSS1P4 = chsSS[0]->p4();
  const LorentzVector& chSS2P4 = chsSS[1]->p4();

  LorentzVector tauP4_trf   = boostToHelicityFrame(tauP4);
  LorentzVector chOSP4_trf  = boostToHelicityFrame(chOSP4);
  LorentzVector chSS1P4_trf = boostToHelicityFrame(chSS1P4);
  LorentzVector chSS2P4_trf = boostToHelicityFrame(chSS2P4);
  reco::Candidate::Vector h = a1pol(chSS1P4_trf, chSS2P4_trf, chOSP4_trf, charge_sum, PolarimetricVectorTau2a1::k3ChargedPi);
  return h.unit();
}
