#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include <TMath.h>
#include <TArrayF.h>
#include <Math/VectorUtil.h>

#include <math.h>

using namespace classic_svFit;

/// global function pointer, needed for Markov Chain integration
const ClassicSVfitIntegrand* ClassicSVfitIntegrand::gSVfitIntegrand = 0;

ClassicSVfitIntegrand::ClassicSVfitIntegrand(int verbosity)
  : beamAxis_(0., 0., 1.),
    invCovMET_(2,2),
#ifdef USE_SVFITTF
    hadTauTF1_(0),
    hadTauTF2_(0),
    useHadTauTF_(false),
    rhoHadTau_(0.),
#endif
    idxLeg1_X_(-1),
    idxLeg1_phi_(-1),
    idxLeg1VisPtShift_(-1),
    idxLeg1_mNuNu_(-1),
    idxLeg2_X_(-1),
    idxLeg2_phi_(-1),
    idxLeg2VisPtShift_(-1),
    idxLeg2_mNuNu_(-1),
    numDimensions_(0),
    addLogM_fixed_(true),
    addLogM_fixed_power_(6.), // CV: best compatibility with "old" SVfitStandalone algorithm
    addLogM_dynamic_(false),
    addLogM_dynamic_formula_(0),
    errorCode_(0),
    histogramAdapter_(0),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<ClassicSVfitIntegrand::ClassicSVfitIntegrand>:" << std::endl;
  }

  // set global function pointer to this
  gSVfitIntegrand = this;
}

ClassicSVfitIntegrand::~ClassicSVfitIntegrand()
{
  if ( verbosity_ ) {
    std::cout << "<ClassicSVfitIntegrand::~ClassicSVfitIntegrand>:" << std::endl;
  }

#ifdef USE_SVFITTF
  delete hadTauTF1_;
  delete hadTauTF2_;
#endif

  delete addLogM_dynamic_formula_;
}

void ClassicSVfitIntegrand::addLogM_fixed(bool value, double power)
{
  addLogM_fixed_ = value;
  addLogM_fixed_power_ = power;
  if ( addLogM_fixed_ && addLogM_dynamic_ ) {
    std::cerr << "Warning: simultaneous use of fixed and dynamic logM terms not supported --> disabling dynamic logM term !!" << std::endl;
    addLogM_dynamic_ = false;
  }
}
void ClassicSVfitIntegrand::addLogM_dynamic(bool value, const std::string& power)
{
  addLogM_dynamic_ = value;
  if ( addLogM_dynamic_ ) {
    if ( power != "" ) {
      TString power_tstring = power.data();
      power_tstring = power_tstring.ReplaceAll("m", "x");
      power_tstring = power_tstring.ReplaceAll("mass", "x");
      std::string formulaName = "ClassicSVfitIntegrand_addLogM_dynamic_formula";
      delete addLogM_dynamic_formula_;
      addLogM_dynamic_formula_ = new TFormula(formulaName.data(), power_tstring.Data());
    } else {
      std::cerr << "Warning: expression = '" << power << "' is invalid --> disabling dynamic logM term !!" << std::endl;
      addLogM_dynamic_ = false;
    }
  }
  if ( addLogM_dynamic_ && addLogM_fixed_ ) {
    std::cerr << "Warning: simultaneous use of fixed and dynamic logM terms not supported --> disabling fixed logM term !!" << std::endl;
    addLogM_fixed_ = false;
  }
}

void ClassicSVfitIntegrand::setDiTauMassConstraint(double diTauMass)
{
  diTauMassConstraint_ = diTauMass;
}

void ClassicSVfitIntegrand::setHistogramAdapter(HistogramAdapter* histogramAdapter)
{
  histogramAdapter_ = histogramAdapter;
}

void ClassicSVfitIntegrand::setIdxLeg1_X(int idx) { idxLeg1_X_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg1_phi(int idx) { idxLeg1_phi_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg1VisPtShift(int idx) { idxLeg1VisPtShift_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg1_mNuNu(int idx) { idxLeg1_mNuNu_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg2_X(int idx) { idxLeg2_X_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg2_phi(int idx) { idxLeg2_phi_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg2VisPtShift(int idx) { idxLeg2VisPtShift_ = idx; }
void ClassicSVfitIntegrand::setIdxLeg2_mNuNu(int idx) { idxLeg2_mNuNu_ = idx; }
void ClassicSVfitIntegrand::setNumDimensions(unsigned numDimensions) { numDimensions_ = numDimensions; }

#ifdef USE_SVFITTF
void ClassicSVfitIntegrand::setHadTauTF(const HadTauTFBase* hadTauTF)
{
  delete hadTauTF1_;
  hadTauTF1_ = hadTauTF->Clone("leg1");
  delete hadTauTF2_;
  hadTauTF2_ = hadTauTF->Clone("leg2");
}

void ClassicSVfitIntegrand::enableHadTauTF()
{
  if ( !(hadTauTF1_ && hadTauTF2_) ) {
    std::cerr << "No tau pT transfer functions defined, call 'setHadTauTF' function first !!" << std::endl;
    assert(0);
  }
  useHadTauTF_ = true;
}
void ClassicSVfitIntegrand::disableHadTauTF()
{
  useHadTauTF_ = false;
}

void ClassicSVfitIntegrand::setRhoHadTau(double rhoHadTau)
{
  rhoHadTau_ = rhoHadTau;
}
#endif


namespace
{
  double norm(const Vector& v)
  {
    return TMath::Sqrt(v.mag2());
  }
}

void
ClassicSVfitIntegrand::setInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET)
{
  if ( verbosity_ ) {
    std::cout << "<ClassicSVfitIntegrand::setInputs>:" << std::endl;
  }

  // reset 'LeptonNumber' and 'MatrixInversion' error codes
  errorCode_ &= (errorCode_ ^ LeptonNumber);
  errorCode_ &= (errorCode_ ^ MatrixInversion);

  if ( measuredTauLeptons.size() != 2 ) {
    std::cerr << "Error: Number of MeasuredTauLeptons is not equal to two !!" << std::endl;
    errorCode_ |= LeptonNumber;
  }
  measuredTauLepton1_ = measuredTauLeptons[0];
  leg1isLep_ = measuredTauLepton1_.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton1_.type() == MeasuredTauLepton::kTauToMuDecay;
  leg1Mass_ = measuredTauLepton1_.mass();
  leg1Mass2_ = square(leg1Mass_);
  Vector eZ1 = normalize(measuredTauLepton1_.p3());
  Vector eY1 = normalize(compCrossProduct(eZ1, beamAxis_));
  Vector eX1 = normalize(compCrossProduct(eY1, eZ1));
  if ( verbosity_ >= 2 ) {
    std::cout << "eX1: theta = " << eX1.theta() << ", phi = " << eX1.phi() << ", norm = " << norm(eX1) << std::endl;
    std::cout << "eY1: theta = " << eY1.theta() << ", phi = " << eY1.phi() << ", norm = " << norm(eY1) << std::endl;
    std::cout << "eZ1: theta = " << eZ1.theta() << ", phi = " << eZ1.phi() << ", norm = " << norm(eZ1) << std::endl;
    std::cout << "(eX1 x eY1 = " << norm(compCrossProduct(eX1, eY1)) << ", eX1 x eZ1 = " << norm(compCrossProduct(eY1, eZ1)) << ", eY1 x eZ1 = " << norm(compCrossProduct(eY1, eZ1)) << ")" << std::endl;
  }
  leg1eX_x_ = eX1.x();
  leg1eX_y_ = eX1.y();
  leg1eX_z_ = eX1.z();
  leg1eY_x_ = eY1.x();
  leg1eY_y_ = eY1.y();
  leg1eY_z_ = eY1.z();
  leg1eZ_x_ = eZ1.x();
  leg1eZ_y_ = eZ1.y();
  leg1eZ_z_ = eZ1.z();

  measuredTauLepton2_ = measuredTauLeptons[1];
  leg2isLep_ = measuredTauLepton2_.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2_.type() == MeasuredTauLepton::kTauToMuDecay;
  leg2Mass_ = measuredTauLepton2_.mass();
  leg2Mass2_ = square(leg2Mass_);
  Vector eZ2 = normalize(measuredTauLepton2_.p3());
  Vector eY2 = normalize(compCrossProduct(eZ2, beamAxis_));
  Vector eX2 = normalize(compCrossProduct(eY2, eZ2));
  if ( verbosity_ >= 2 ) {
    std::cout << "eX2: theta = " << eX2.theta() << ", phi = " << eX2.phi() << ", norm = " << norm(eX2) << std::endl;
    std::cout << "eY2: theta = " << eY2.theta() << ", phi = " << eY2.phi() << ", norm = " << norm(eY2) << std::endl;
    std::cout << "eZ2: theta = " << eZ2.theta() << ", phi = " << eZ2.phi() << ", norm = " << norm(eZ2) << std::endl;
    //std::cout << "(eX2 x eY2 = " << norm(compCrossProduct(eX2, eY2)) << ", eX2 x eZ2 = " << norm(compCrossProduct(eY2, eZ2)) << ", eY2 x eZ2 = " << norm(compCrossProduct(eY2, eZ2)) << ")" << std::endl;
  }
  leg2eX_x_ = eX2.x();
  leg2eX_y_ = eX2.y();
  leg2eX_z_ = eX2.z();
  leg2eY_x_ = eY2.x();
  leg2eY_y_ = eY2.y();
  leg2eY_z_ = eY2.z();
  leg2eZ_x_ = eZ2.x();
  leg2eZ_y_ = eZ2.y();
  leg2eZ_z_ = eZ2.z();

  mVis_measured_ = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).mass();
  if ( verbosity_ >= 2 ) {
    std::cout << "mVis = " << mVis_measured_ << std::endl;
  }
  mVis2_measured_ = square(mVis_measured_);

  measuredMETx_ = measuredMETx;
  measuredMETy_ = measuredMETy;

  // determine transfer matrix for MET
  invCovMET_ = covMET;
  double covDet = invCovMET_.Determinant();
  const_MET_ = 0.;
  if ( covDet != 0 ) {
    invCovMET_.Invert();
    invCovMETxx_ = invCovMET_(0,0);
    invCovMETxy_ = invCovMET_(0,1);
    invCovMETyx_ = invCovMET_(1,0);
    invCovMETyy_ = invCovMET_(1,1);
    const_MET_ = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));
  } else{
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!" << std::endl;
    errorCode_ |= MatrixInversion;
  }

#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) {
    if ( measuredTauLepton1_.type() == MeasuredTauLepton::kTauToHadDecay ) {
      assert(hadTauTF1_);
      hadTauTF1_->setDecayMode(measuredTauLepton1_.decayMode());
    }
    if ( measuredTauLepton2_.type() == MeasuredTauLepton::kTauToHadDecay ) {
      assert(hadTauTF2_);
      hadTauTF2_->setDecayMode(measuredTauLepton2_.decayMode());
    }
  }
#endif
}

double
ClassicSVfitIntegrand::Eval(const double* x) const
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand::Eval(const double*)>:" << std::endl;
    std::cout << " x = { ";
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << x[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ != 0 ) {
    return 0.;
  }

  double visPtShift1 = ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) ? (1./x[idxLeg1VisPtShift_]) : 1.;
  double visPtShift2 = ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) ? (1./x[idxLeg2VisPtShift_]) : 1.;
  if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 ) return 0.;

  // compute four-vector of visible decay products for first tau
  double vis1Px = visPtShift1*measuredTauLepton1_.px();
  double vis1Py = visPtShift1*measuredTauLepton1_.py();
  double vis1Pz = visPtShift1*measuredTauLepton1_.pz();
  double vis1En = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  //std::cout << "vis1: En = " << vis1En << ", Pt = " << TMath::Sqrt(vis1Px*vis1Px + vis1Py*vis1Py) << std::endl;
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  double vis1P = vis1P4.P();

  // compute four-vector of visible decay products for second tau
  double vis2Px = visPtShift2*measuredTauLepton2_.px();
  double vis2Py = visPtShift2*measuredTauLepton2_.py();
  double vis2Pz = visPtShift2*measuredTauLepton2_.pz();
  double vis2En = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  //std::cout << "vis2: En = " << vis2En << ", Pt = " << TMath::Sqrt(vis2Px*vis2Px + vis2Py*vis2Py) << std::endl;
  LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  double vis2P = vis2P4.P();

  // compute visible energy fractions for both taus
  assert(idxLeg1_X_ != -1);
  double x1_dash = x[idxLeg1_X_];
  double x1 = x1_dash/visPtShift1;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  double x2_dash = 0.0;
  if (idxLeg2_X_ != -1) {
    x2_dash = x[idxLeg2_X_];
  }
  else {
    x2_dash = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).M2()/(diTauMassConstraint_ * diTauMassConstraint_)/x1_dash;
  }
  double x2 = x2_dash/visPtShift2;
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double nu1En = vis1En*(1. - x1)/x1;
  double nu1Mass = ( idxLeg1_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg1_mNuNu_]) : 0.;
  double nu1P = TMath::Sqrt(TMath::Max(0., nu1En*nu1En - square(nu1Mass)));
  assert(idxLeg1_phi_ != -1);
  double phiNu1 = x[idxLeg1_phi_];
  double cosThetaNu1 = compCosThetaNuNu(vis1En, vis1P, leg1Mass2_, nu1En, nu1P, square(nu1Mass));
  if ( !(cosThetaNu1 >= -1. && cosThetaNu1 <= +1.) ) return 0.;
  double thetaNu1 = TMath::ACos(cosThetaNu1);
  double nu1Px_local = nu1P*TMath::Cos(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Py_local = nu1P*TMath::Sin(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Pz_local = nu1P*TMath::Cos(thetaNu1);
  double nu1Px = nu1Px_local*leg1eX_x_ + nu1Py_local*leg1eY_x_ + nu1Pz_local*leg1eZ_x_;
  double nu1Py = nu1Px_local*leg1eX_y_ + nu1Py_local*leg1eY_y_ + nu1Pz_local*leg1eZ_y_;
  double nu1Pz = nu1Px_local*leg1eX_z_ + nu1Py_local*leg1eY_z_ + nu1Pz_local*leg1eZ_z_;
  //std::cout << "nu1: En = " << nu1En << ", Pt = " << TMath::Sqrt(nu1Px*nu1Px + nu1Py*nu1Py) << std::endl;
  LorentzVector nu1P4(nu1Px, nu1Py, nu1Pz, nu1En);

  double tau1En = vis1En + nu1En;
  double tau1Px = vis1P4.px() + nu1Px;
  double tau1Py = vis1P4.py() + nu1Py;
  double tau1Pz = vis1P4.pz() + nu1Pz;
  //std::cout << "tau1: En = " << tau1En << ", Pt = " << TMath::Sqrt(tau1Px*tau1Px + tau1Py*tau1Py) << std::endl;
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);

  // compute neutrino and tau lepton four-vector for second tau
  double nu2En = vis2En*(1. - x2)/x2;
  double nu2Mass = ( idxLeg2_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg2_mNuNu_]) : 0.;
  double nu2P = TMath::Sqrt(TMath::Max(0., nu2En*nu2En - square(nu2Mass)));
  assert(idxLeg2_phi_ != -2);
  double phiNu2 = x[idxLeg2_phi_];
  double cosThetaNu2 = compCosThetaNuNu(vis2En, vis2P, leg2Mass2_, nu2En, nu2P, square(nu2Mass));
  if ( !(cosThetaNu2 >= -1. && cosThetaNu2 <= +1.) ) return 0.;
  double thetaNu2 = TMath::ACos(cosThetaNu2);
  double nu2Px_local = nu2P*TMath::Cos(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Py_local = nu2P*TMath::Sin(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Pz_local = nu2P*TMath::Cos(thetaNu2);
  double nu2Px = nu2Px_local*leg2eX_x_ + nu2Py_local*leg2eY_x_ + nu2Pz_local*leg2eZ_x_;
  double nu2Py = nu2Px_local*leg2eX_y_ + nu2Py_local*leg2eY_y_ + nu2Pz_local*leg2eZ_y_;
  double nu2Pz = nu2Px_local*leg2eX_z_ + nu2Py_local*leg2eY_z_ + nu2Pz_local*leg2eZ_z_;
  //std::cout << "nu2: En = " << nu2En << ", Pt = " << TMath::Sqrt(nu2Px*nu2Px + nu2Py*nu2Py) << std::endl;
  LorentzVector nu2P4(nu2Px, nu2Py, nu2Pz, nu2En);

  double tau2En = vis2En + nu2En;
  double tau2Px = vis2P4.px() + nu2Px;
  double tau2Py = vis2P4.py() + nu2Py;
  double tau2Pz = vis2P4.pz() + nu2Pz;
  //std::cout << "tau2: En = " << tau2En << ", Pt = " << TMath::Sqrt(tau2Px*tau2Px + tau2Py*tau2Py) << std::endl;
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);

  // evaluate transfer function for MET/hadronic recoil
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double residualX = measuredMETx_ - sumNuPx;
  double residualY = measuredMETy_ - sumNuPy;
#ifdef USE_SVFITTF
  if ( rhoHadTau_ != 0. ) {
    residualX += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.px() + (visPtShift2 - 1.)*measuredTauLepton2_.px()));
    residualY += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.py() + (visPtShift2 - 1.)*measuredTauLepton2_.py()));
  }
#endif
  double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
  double prob_TF_met = const_MET_*TMath::Exp(-0.5*pull2);
  if ( verbosity_ >= 2 ) {
    std::cout << "TF(met): recPx = " << measuredMETx_ << ", recPy = " << measuredMETy_ << ", genPx = " << sumNuPx << ", genPy = " << sumNuPy << " --> prob = " << prob_TF_met << std::endl;
    std::cout << "leg1: En = " << vis1P4.energy() << ", Px = " << vis1P4.px() << ", Py = " << vis1P4.py() << ", Pz = " << vis1P4.pz() << ";"
              << " Pt = " << vis1P4.pt() << ", eta = " << vis1P4.eta() << ", phi = " << vis1P4.phi() << ", mass = " << vis1P4.mass()
              << " (x = " << x1 << ")" << std::endl;
    std::cout << "tau1: En = " << tau1P4.energy() << ", Px = " << tau1P4.px() << ", Py = " << tau1P4.py() << ", Pz = " << tau1P4.pz() << ";"
              << " Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << std::endl;
    std::cout << "nu1: En = " << nu1P4.energy() << ", Px = " << nu1P4.px() << ", Py = " << nu1P4.py() << ", Pz = " << nu1P4.pz() << ";"
              << " Pt = " << nu1P4.pt() << ", eta = " << nu1P4.eta() << ", phi = " << nu1P4.phi() << ", mass = " << nu1P4.mass() << std::endl;
    //double angle1 = compAngle(vis1P4, nu1P4);
    //std::cout << "angle(vis1, nu1) = " << angle1 << std::endl;
    //double phiInvis1 = compPhiInvis(vis1P4, nu1P4);
    //std::cout << "phiInvis1 = " << phiInvis1 << std::endl;
    std::cout << "leg2: En = " << vis2P4.energy() << ", Px = " << vis2P4.px() << ", Py = " << vis2P4.py() << ", Pz = " << vis2P4.pz() << ";"
              << " Pt = " << vis2P4.pt() << ", eta = " << vis2P4.eta() << ", phi = " << vis2P4.phi() << ", mass = " << vis2P4.mass()
              << " (x = " << x2 << ")" << std::endl;
    std::cout << "tau2: En = " << tau2P4.energy() << ", Px = " << tau2P4.px() << ", Py = " << tau2P4.py() << ", Pz = " << tau2P4.pz() << ";"
              << " Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << std::endl;
    std::cout << "nu2: En = " << nu2P4.energy() << ", Px = " << nu2P4.px() << ", Py = " << nu2P4.py() << ", Pz = " << nu2P4.pz() << ";"
              << " Pt = " << nu2P4.pt() << ", eta = " << nu2P4.eta() << ", phi = " << nu2P4.phi() << ", mass = " << nu2P4.mass() << std::endl;
    //double angle2 = compAngle(vis2P4, nu2P4);
    //std::cout << "angle(vis2, nu2) = " << angle2 << std::endl;
    //double phiInvis2 = compPhiInvis(vis2P4, nu2P4);
    //std::cout << "phiInvis2 = " << phiInvis2 << std::endl;
  }
  double prob_TF = prob_TF_met;

#ifdef USE_SVFITTF
  // evaluate transfer functions for tau energy reconstruction
  if ( useHadTauTF_ && idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) {
    double prob_TF_leg1 = (*hadTauTF1_)(measuredTauLepton1_.pt(), vis1P4.pt(), vis1P4.eta());
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg1): recPt = " << measuredTauLepton1_.pt() << ", genPt = " << vis1P4.pt() << ", genEta = " << vis1P4.eta() << " --> prob = " << prob_TF_leg1 << std::endl;
    }
    prob_TF *= prob_TF_leg1;
  }
  if ( useHadTauTF_ && idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) {
    double prob_TF_leg2 = (*hadTauTF2_)(measuredTauLepton2_.pt(), vis2P4.pt(), vis2P4.eta());
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg2): recPt = " << measuredTauLepton2_.pt() << ", genPt = " << vis2P4.pt() << ", genEta = " << vis2P4.eta() << " --> prob = " << prob_TF_leg2 << std::endl;
    }
    prob_TF *= prob_TF_leg2;
  }
#endif

  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m//FIX ME store this
  const double constFactor = 2.*conversionFactor/eigth(2.*TMath::Pi()); //FIXME store this
  double prob_PS_and_tauDecay = constFactor;
  double prob_tauDecay_leg1 = 0.;
  if ( leg1isLep_ ) {
    prob_tauDecay_leg1 = compPSfactor_tauToLepDecay(x1, vis1P4.E(), vis1P4.P(), leg1Mass_, nu1P4.E(), nu1P4.P(), nu1Mass);
  } else {
    prob_tauDecay_leg1 = compPSfactor_tauToHadDecay(x1, vis1P4.E(), vis1P4.P(), leg1Mass_, nu1P4.E(), nu1P4.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg1;
  double prob_tauDecay_leg2 = 0.;
  if ( leg2isLep_ ) {
    prob_tauDecay_leg2 = compPSfactor_tauToLepDecay(x2, vis2P4.E(), vis2P4.P(), leg2Mass_, nu2P4.E(), nu2P4.P(), nu2Mass);
  } else {
    prob_tauDecay_leg2 = compPSfactor_tauToHadDecay(x2, vis2P4.E(), vis2P4.P(), leg2Mass_, nu2P4.E(), nu2P4.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg2;
  // CV: multiply matrix element by factor (Pi/(mTau GammaTau))^2 from Luca's write-up
  prob_PS_and_tauDecay *= square(TMath::Pi()/(tauLeptonMass*GammaTau));

  double mTauTau = (tau1P4 + tau2P4).mass();
  double prob_logM = 1.;
  if ( addLogM_fixed_ ) {
    prob_logM = 1./TMath::Power(TMath::Max(1., mTauTau), addLogM_fixed_power_);
  }
  if ( addLogM_dynamic_ ) {
    double addLogM_power = addLogM_dynamic_formula_->Eval(mTauTau);
    prob_logM = 1./TMath::Power(TMath::Max(1., mTauTau), TMath::Max(0., addLogM_power));
  }

  double jacobiFactor = 1./(visPtShift1*visPtShift2); // product of derrivatives dx1/dx1' and dx2/dx2' for parametrization of x1, x2 by x1', x2'
  if (diTauMassConstraint_ > 0.0) {
    jacobiFactor *= (2.0*x2*diTauMassConstraint_);
  }
  double prob = prob_PS_and_tauDecay*prob_TF*prob_logM*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: PS+decay = " << prob_PS_and_tauDecay << ","
              << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  }
  if ( TMath::IsNaN(prob) ) {
    std::cerr << "Warning: prob = " << prob << " (PS+decay = " << prob_PS_and_tauDecay << ","
              << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor << ") --> setting prob = 0 !!" << std::endl;
    prob = 0.;
  }

  // CV: fill histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  histogramAdapter_->setTau1P4(tau1P4);
  histogramAdapter_->setTau2P4(tau2P4);

  return prob;
}
