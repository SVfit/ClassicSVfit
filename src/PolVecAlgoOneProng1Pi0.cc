#include "TauAnalysis/ClassicSVfit/interface/PolVecAlgoOneProng1Pi0.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // fixNuMass(), gamma_va, tauLeptonMass

#include <assert.h>                                               // assert()

namespace classic_svFit
{

PolVecAlgoOneProng1Pi0::PolVecAlgoOneProng1Pi0()
{}

PolVecAlgoOneProng1Pi0::~PolVecAlgoOneProng1Pi0()
{}

classic_svFit::Vector
PolVecAlgoOneProng1Pi0::operator()(const MeasuredTauLepton& measuredTauLepton, const FittedTauLepton& fittedTauLepton, int tau, 
                                   const BoostToHelicityFrame& boostToHelicityFrame) const
{
  const std::vector<MeasuredHadTauDecayProduct>& daughters = measuredTauLepton.hadTauDecayProducts();
  const MeasuredHadTauDecayProduct* ch = nullptr;
  const MeasuredHadTauDecayProduct* pi0 = nullptr;
  for ( const MeasuredHadTauDecayProduct& daughter : daughters )
  {
    if ( daughter.charge() != 0 )
    {
      ch = &daughter;
    }
    else
    {
      pi0 = &daughter;
    }
  }
  if ( !ch )
  {
    std::cerr << "ERROR: Failed to find charged pion !!" << std::endl;
    assert(0);
  }
  if ( !pi0 )
  {
    std::cerr << "ERROR: Failed to find neutral pion !!" << std::endl;
    assert(0);
  }

  // CV: notation of four-vectors chosen according to Section 3.4 of the paper
  //       Comput.Phys.Commun. 64 (1990) 275
  LorentzVector q1 = boostToHelicityFrame(ch->p4(), tau);
  LorentzVector q2 = boostToHelicityFrame(pi0->p4(), tau);

  // CV: the neutrino four-vector is computed by taking the difference 
  //     between the tau four-vector and the four-vector of the visible tau decay products
  //     yields mass values that a few GeV off, presumably due to rounding errors.
  //     Adjust the energy component of the neutrino four-vector such that its mass value equals zero,
  //     while keeping the Px, Py, Pz momentum components fixed
  const LorentzVector& tauP4 = fittedTauLepton.tauP4();
  const LorentzVector& visP4 = fittedTauLepton.visP4();
  LorentzVector nuP4 = fixNuMass(tauP4 - visP4);
  LorentzVector N = fixNuMass(boostToHelicityFrame(nuP4, tau));
  assert(nuP4.energy() >= 0. && N.energy() >= 0.);

  LorentzVector P = boostToHelicityFrame(tauP4, tau);

  LorentzVector q = q1 - q2;
  double omega = 2.*(q.Dot(N))*(q.Dot(P)) - q.mass2()*(N.Dot(P));
  // CV: sign of terms proportional to gammaVA differs for tau+ and tau-,
  //     cf. text following Eq. (3.16) in Comput.Phys.Commun. 64 (1991) 275
  double sign = 0.;
  if      ( ch->charge() > 0. ) sign = -1.;
  else if ( ch->charge() < 0. ) sign = +1.;
  else assert(0);
  // CV: term 2.*|f2|^2 appears in expression for h as well as in expression for omega
  //     and drops out
  Vector h = -(sign*gamma_va*tauLeptonMass/omega)*(2.*(q.Dot(N))*q.Vect() - q.mass2()*N.Vect());
  return h.unit();
}

}
