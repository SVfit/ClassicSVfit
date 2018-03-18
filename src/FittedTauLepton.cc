#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"

#include <TMath.h>

using namespace classic_svFit;

FittedTauLepton::FittedTauLepton(int iTau, int verbosity)
  : iTau_(iTau)
  , measuredTauLepton_mass_(-1.)
  , measuredTauLepton_mass2_(-1.)
  , x_(-1.)
  , phiNu_(-1.)
  , nuMass_(-1.)
  , errorCode_(None)
  , verbosity_(verbosity)
{}

FittedTauLepton::~FittedTauLepton()
{}

namespace
{
  double norm(const Vector& v)
  {
    return TMath::Sqrt(v.mag2());
  }
}

void FittedTauLepton::setMeasuredTauLepton(const MeasuredTauLepton& measuredTauLepton)
{
  measuredTauLepton_ = measuredTauLepton;
  measuredTauLepton_mass_ = measuredTauLepton.mass();
  measuredTauLepton_mass2_ = square(measuredTauLepton_mass_);

  Vector beamAxis(0., 0., 1.);
  Vector eZ = normalize(measuredTauLepton_.p3());
  Vector eY = normalize(compCrossProduct(eZ, beamAxis));
  Vector eX = normalize(compCrossProduct(eY, eZ));
  if ( verbosity_ >= 2 ) {
    std::cout << "eX" << (iTau_ + 1) << ": theta = " << eX.theta() << ", phi = " << eX.phi() << ", norm = " << norm(eX) << std::endl;
    std::cout << "eY" << (iTau_ + 1) << ": theta = " << eY.theta() << ", phi = " << eY.phi() << ", norm = " << norm(eY) << std::endl;
    std::cout << "eZ" << (iTau_ + 1) << ": theta = " << eZ.theta() << ", phi = " << eZ.phi() << ", norm = " << norm(eZ) << std::endl;
    std::cout << "(eX" << (iTau_ + 1) << " x eY" << (iTau_ + 1) << " = " << norm(compCrossProduct(eX, eY)) << "," 
	      << " eX" << (iTau_ + 1) << " x eZ" << (iTau_ + 1) << " = " << norm(compCrossProduct(eY, eZ)) << "," 
	      << " eY" << (iTau_ + 1) << " x eZ" << (iTau_ + 1) << " = " << norm(compCrossProduct(eY, eZ)) << ")" << std::endl;
  }
  eX_x_ = eX.x();
  eX_y_ = eX.y();
  eX_z_ = eX.z();
  eY_x_ = eY.x();
  eY_y_ = eY.y();
  eY_z_ = eY.z();
  eZ_x_ = eZ.x();
  eZ_y_ = eZ.y();
  eZ_z_ = eZ.z();
}

const MeasuredTauLepton& FittedTauLepton::getMeasuredTauLepton() const
{
  return measuredTauLepton_;
}

void FittedTauLepton::updateVisMomentum(double visPtShift)
{
  // compute four-vector of visible decay products
  double visPx = visPtShift*measuredTauLepton_.px();
  double visPy = visPtShift*measuredTauLepton_.py();
  double visPz = visPtShift*measuredTauLepton_.pz();
  double visEn = TMath::Sqrt(square(visPx) + square(visPy) + square(visPz) + measuredTauLepton_mass2_);
  //std::cout << "vis: En = " << visEn << ", Pt = " << TMath::Sqrt(square(visPx) + square(visPy)) << std::endl;
  visP4_.SetPxPyPzE(visPx, visPy, visPz, visEn);

  // set tau lepton four-vector to four-vector of visible decay products and neutrino four-vector to zero,
  // in case of electrons or muons directly originating from LFV Higgs boson decay
  if ( measuredTauLepton_.type() == MeasuredTauLepton::kPrompt ) {
    nuP4_.SetPxPyPzE(0., 0., 0., 0.);
    tauP4_ = visP4_;
  }
}

void FittedTauLepton::updateTauMomentum(double x, double phiNu, double nuMass)
{
  x_ = x;
  phiNu_ = phiNu;
  nuMass_ = nuMass;

  errorCode_ = None;

  // compute neutrino and tau lepton four-vector 
  double nuEn = visP4_.E()*(1. - x_)/x_;
  double nuMass2 = square(nuMass_);
  double nuP = TMath::Sqrt(TMath::Max(0., square(nuEn) - nuMass2));
  double cosThetaNu = compCosThetaNuNu(visP4_.E(), visP4_.P(), measuredTauLepton_mass2_, nuEn, nuP, nuMass2);
  if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) {
    errorCode_ |= TauDecayParameters;
    return;
  }

  double cosPhiNu, sinPhiNu;
  sincos(phiNu_, &sinPhiNu, &cosPhiNu);
  double thetaNu = TMath::ACos(cosThetaNu);
  double sinThetaNu = TMath::Sin(thetaNu);

  double nuPx_local = nuP*cosPhiNu*sinThetaNu;
  double nuPy_local = nuP*sinPhiNu*sinThetaNu;
  double nuPz_local = nuP*cosThetaNu;
  double nuPx = nuPx_local*eX_x_ + nuPy_local*eY_x_ + nuPz_local*eZ_x_;
  double nuPy = nuPx_local*eX_y_ + nuPy_local*eY_y_ + nuPz_local*eZ_y_;
  double nuPz = nuPx_local*eX_z_ + nuPy_local*eY_z_ + nuPz_local*eZ_z_;
  //std::cout << "nu1: En = " << nuEn << ", Pt = " << TMath::Sqrt(square(nuPx) + square(nuPy)) << std::endl;
  nuP4_.SetPxPyPzE(nuPx, nuPy, nuPz, nuEn);

  double tauEn = visP4_.E() + nuEn;
  double tauPx = visP4_.px() + nuPx;
  double tauPy = visP4_.py() + nuPy;
  double tauPz = visP4_.pz() + nuPz;
  //std::cout << "tau: En = " << tauEn << ", Pt = " << TMath::Sqrt(square(tauPx) + square(tauPy)) << std::endl;
  tauP4_.SetPxPyPzE(tauPx, tauPy, tauPz, tauEn);
}

const LorentzVector& FittedTauLepton::visP4() const
{
  return visP4_;
}

const LorentzVector& FittedTauLepton::nuP4() const
{
  return nuP4_;
}

const LorentzVector& FittedTauLepton::tauP4() const
{
  return tauP4_;
}

double FittedTauLepton::x() const
{
  return x_;
}

double FittedTauLepton::phiNu() const
{
  return phiNu_;
}

double FittedTauLepton::nuMass() const
{
  return nuMass_;
}

int FittedTauLepton::errorCode() const
{
  return errorCode_;
}
