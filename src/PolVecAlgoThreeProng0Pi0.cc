#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoThreeProng0Pi0.h"

#include <assert.h> // assert()

namespace classic_svFit
{

PolVecAlgoThreeProng0Pi0::PolVecAlgoThreeProng0Pi0()
{}

PolVecAlgoThreeProng0Pi0::~PolVecAlgoThreeProng0Pi0()
{}

}

namespace
{
  bool
  isHigherPt(const classic_svFit::MeasuredHadTauDecayProduct* particle1, const classic_svFit::MeasuredHadTauDecayProduct* particle2)
  {
    return particle1->pt() > particle2->pt();
  }
}

namespace classic_svFit
{

Vector
PolVecAlgoThreeProng0Pi0::operator()(const MeasuredTauLepton& measuredTauLepton, int tau,
                                     const BoostToHelicityFrame& boostToHelicityFrame) const
{
  const std::vector<MeasuredHadTauDecayProduct>& daughters = measuredTauLepton.hadTauDecayProducts();
  std::vector<const MeasuredHadTauDecayProduct*> chsPlus;
  std::vector<const MeasuredHadTauDecayProduct*> chsMinus;
  int charge_sum = 0;
  for ( const MeasuredHadTauDecayProduct& daughter : daughters )
  {
    if      ( daughter.charge() > 0 ) chsPlus.push_back(&daughter);
    else if ( daughter.charge() < 0 ) chsMinus.push_back(&daughter);
    charge_sum += daughter.charge();
  }
  if ( !((chsPlus.size() + chsMinus.size()) == 3 && std::abs(charge_sum) == 1) )
  {
    std::cerr << "ERROR: Failed to find three charged pions with |sum(charge)| = 1 !!" << std::endl;
    assert(0);
  }

  std::vector<const MeasuredHadTauDecayProduct*> chsOS;
  std::vector<const MeasuredHadTauDecayProduct*> chsSS;
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

  LorentzVector chOSP4_trf  = boostToHelicityFrame(chOSP4, tau);
  LorentzVector chSS1P4_trf = boostToHelicityFrame(chSS1P4, tau);
  LorentzVector chSS2P4_trf = boostToHelicityFrame(chSS2P4, tau);
  Vector h = a1pol_(chSS1P4_trf, chSS2P4_trf, chOSP4_trf, charge_sum, PolarimetricVectorTau2a1::k3ChargedPi);
  return h.unit();
}

}
