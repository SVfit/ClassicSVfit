#include "TauAnalysis/ClassicSVfit/interface/BoostToHelicityFrame.h"

#include <vector> // std::vector<>

namespace
{
  classic_svFit::LorentzVector
  getP4_rf(const classic_svFit::LorentzVector& p4, const ROOT::Math::Boost& boost)
  {
    // CV: boost given four-vector to restframe
    classic_svFit::LorentzVector p4_rf = boost(p4);
    return p4_rf;
  }

  classic_svFit::Vector
  get_k(const classic_svFit::LorentzVector& p4, const ROOT::Math::Boost& boost_ttrf)
  {
    classic_svFit::LorentzVector p4_ttrf = getP4_rf(p4, boost_ttrf);
    classic_svFit::Vector k = p4_ttrf.Vect().unit();
    return k;
  }

  classic_svFit::Vector
  get_h(const classic_svFit::LorentzVector& beamP4, const ROOT::Math::Boost& boost_ttrf)
  {
    classic_svFit::LorentzVector beamP4_ttrf = getP4_rf(beamP4, boost_ttrf);
    classic_svFit::Vector h = beamP4_ttrf.Vect().unit();
    return h;
  }

  void
  get_localCoordinateSystem(const classic_svFit::LorentzVector& tauP4, const ROOT::Math::Boost& boost_ttrf,
                            const classic_svFit::LorentzVector& beamP4, 
                            classic_svFit::Vector& r, classic_svFit::Vector& n, classic_svFit::Vector& k)
  {
    k = get_k(tauP4, boost_ttrf);
    classic_svFit::Vector h = get_h(beamP4, boost_ttrf);
    r = classic_svFit::get_r(k, h);
    n = classic_svFit::get_n(k, r);
  }

  classic_svFit::LorentzVector
  getP4_hf(const classic_svFit::LorentzVector& p4, const classic_svFit::Vector& r, const classic_svFit::Vector& n, const classic_svFit::Vector& k)
  {
    // CV: rotate given four-vector to helicity frame
    classic_svFit::Vector p3 = p4.Vect();
    double Pr = p3.Dot(r);
    double Pn = p3.Dot(n);
    double Pk = p3.Dot(k);
    classic_svFit::LorentzVector p4_hf(Pr, Pn, Pk, p4.energy());
    return p4_hf;
  }
}

namespace classic_svFit
{

BoostToHelicityFrame::BoostToHelicityFrame()
{
  beamP4_ = get_beamP4();
}

BoostToHelicityFrame::~BoostToHelicityFrame()
{}

void
BoostToHelicityFrame::setFittedTauLeptons(const FittedTauLepton& fittedTauLepton1, const FittedTauLepton& fittedTauLepton2)
{
  std::vector<const FittedTauLepton*> fittedTauLeptons;
  fittedTauLeptons.push_back(&fittedTauLepton1);
  fittedTauLeptons.push_back(&fittedTauLepton2);
  const FittedTauLepton* fittedTauPlus = nullptr;
  const FittedTauLepton* fittedTauMinus = nullptr;
  for ( const FittedTauLepton* fittedTauLepton : fittedTauLeptons )
  {
    if ( fittedTauLepton->getMeasuredTauLepton().charge() > 0 )
    {
      fittedTauPlus = fittedTauLepton;
    }
    else if ( fittedTauLepton->getMeasuredTauLepton().charge() < 0 )
    {
      fittedTauMinus = fittedTauLepton;
    }
  }
  assert(fittedTauPlus && fittedTauMinus);
  LorentzVector tauPlusP4 = fittedTauPlus->tauP4();
  LorentzVector tauMinusP4 = fittedTauMinus->tauP4();
  LorentzVector diTauP4 = tauPlusP4 + tauMinusP4;

  boost_ttrf_ = ROOT::Math::Boost(diTauP4.BoostToCM());
  ::get_localCoordinateSystem(tauMinusP4, boost_ttrf_, beamP4_, r_, n_, k_);

  LorentzVector tauPlusP4_ttrf = getP4_rf(tauPlusP4, boost_ttrf_);
  LorentzVector tauPlusP4_hf = getP4_hf(tauPlusP4_ttrf, r_, n_, k_);
  boost_tprf_ = ROOT::Math::Boost(tauPlusP4_hf.BoostToCM());

  LorentzVector tauMinusP4_ttrf = getP4_rf(tauMinusP4, boost_ttrf_);
  LorentzVector tauMinusP4_hf = getP4_hf(tauMinusP4_ttrf, r_, n_, k_);
  boost_tmrf_ = ROOT::Math::Boost(tauMinusP4_hf.BoostToCM());
}

LorentzVector
BoostToHelicityFrame::operator()(const LorentzVector& p4, int tau) const
{
  // CV: boost given four-vector to restframe of tau pair,
  //     rotate to helicity frame,
  //     and finally boost to tau restframe
  LorentzVector p4_ttrf = getP4_rf(p4, boost_ttrf_);
  LorentzVector p4_hf = getP4_hf(p4_ttrf, r_, n_, k_);
  const ROOT::Math::Boost* boost_trf = nullptr;
  if      ( tau == kTauPlus  ) boost_trf = &boost_tprf_;
  else if ( tau == kTauMinus ) boost_trf = &boost_tmrf_;
  assert(boost_trf);
  LorentzVector p4_trf = getP4_rf(p4_hf, *boost_trf);
  return p4_trf;
}

}
