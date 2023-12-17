#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoOneProng0Pi0.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // gamma_va, tauLeptonMass2, tauLeptonMass3

#include <assert.h>                                               // assert()

namespace classic_svFit
{

PolVecAlgoOneProng0Pi0::PolVecAlgoOneProng0Pi0()
{}

PolVecAlgoOneProng0Pi0::~PolVecAlgoOneProng0Pi0()
{}

Vector
PolVecAlgoOneProng0Pi0::operator()(const MeasuredTauLepton& measuredTauLepton, int tau,
                                   const BoostToHelicityFrame& boostToHelicityFrame) const
{
  const std::vector<MeasuredHadTauDecayProduct>& daughters = measuredTauLepton.hadTauDecayProducts();
  const MeasuredHadTauDecayProduct* ch = nullptr;
  for ( const MeasuredHadTauDecayProduct& daughter : daughters )
  {
    if ( daughter.charge() != 0 )
    {
      ch = &daughter;
    } 
  }
  if ( !ch )
  {
    std::cerr << "ERROR: Failed to find charged pion !!" << std::endl;
    assert(0);
  }

  // CV: notation of four-vectors chosen according to Section 3.3 of the paper
  //       Comput.Phys.Commun. 64 (1990) 275
  const LorentzVector& chP4 = ch->p4();
  LorentzVector Q = boostToHelicityFrame(chP4, tau);

  const double f1 = 0.1284;
  double omega = (tauLeptonMass2 - chP4.mass2())*tauLeptonMass2;
  // CV: sign of terms proportional to gammaVA differs for tau+ and tau-,
  //     cf. text following Eq. (3.16) in Comput.Phys.Commun. 64 (1991) 275
  double sign = 0.;
  if      ( ch->charge() > 0. ) sign = -1.;
  else if ( ch->charge() < 0. ) sign = +1.;
  else assert(0);

  Vector h = -(2.*sign*gamma_va*square(f1)*tauLeptonMass3/omega)*Q.Vect();
  return h.unit();
}

}
