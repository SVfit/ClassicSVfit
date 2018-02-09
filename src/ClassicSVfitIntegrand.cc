#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"

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
#ifdef USE_SVFITTF
    hadTauTF1_(0),
    hadTauTF2_(0),
    useHadTauTF_(false),
    rhoHadTau_(0.),
#endif
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

  vis1En_ = -1;
  vis2En_ = -1;

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


void ClassicSVfitIntegrand::setVerbosity(int aVerbosity){ verbosity_ = aVerbosity;}

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

void ClassicSVfitIntegrand::setLegIntegrationParams(unsigned int iLeg,
                                                const classic_svFit::integrationParameters & aParams)
                                                { legIntegrationParams_[iLeg] = aParams;}
void ClassicSVfitIntegrand::setNumDimensions(unsigned numDimensions) { numDimensions_ = numDimensions; }

void ClassicSVfitIntegrand::setIntegrationRanges(const double* xl, const double* xu){
 for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
  }
}

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
ClassicSVfitIntegrand::setLeptonInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  if ( verbosity_>=2 ) {
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

  computeVisMom(1.0, 1.0);

  phaseSpaceComponentCache_ = 0;

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

void ClassicSVfitIntegrand::addMETEstimate(double measuredMETx, double measuredMETy, const TMatrixD& covMET){

  measuredMETx_.push_back(measuredMETx);
  measuredMETy_.push_back(measuredMETy);
  covMET_.push_back(covMET);

}

int ClassicSVfitIntegrand::getMETComponentsSize() const {return measuredMETx_.size();}

void ClassicSVfitIntegrand::clearMET(){
  measuredMETx_.clear();
  measuredMETy_.clear();
  covMET_.clear();
}

void ClassicSVfitIntegrand::computeVisMom(const double & visPtShift1, const double & visPtShift2){

// compute four-vector of visible decay products for first tau
  double vis1Px = visPtShift1*measuredTauLepton1_.px();
  double vis1Py = visPtShift1*measuredTauLepton1_.py();
  double vis1Pz = visPtShift1*measuredTauLepton1_.pz();
  vis1En_ = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  //std::cout << "vis1: En = " << vis1En << ", Pt = " << TMath::Sqrt(vis1Px*vis1Px + vis1Py*vis1Py) << std::endl;
  vis1P4_.SetPxPyPzE(vis1Px, vis1Py, vis1Pz, vis1En_);
  vis1P_ = vis1P4_.P();

  // compute four-vector of visible decay products for second tau
  double vis2Px = visPtShift2*measuredTauLepton2_.px();
  double vis2Py = visPtShift2*measuredTauLepton2_.py();
  double vis2Pz = visPtShift2*measuredTauLepton2_.pz();
  vis2En_ = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  //std::cout << "vis2: En = " << vis2En << ", Pt = " << TMath::Sqrt(vis2Px*vis2Px + vis2Py*vis2Py) << std::endl;
  vis2P4_.SetPxPyPzE(vis2Px, vis2Py, vis2Pz, vis2En_);
  vis2P_ = vis2P4_.P();
}

void ClassicSVfitIntegrand::rescaleX(const double* q) const
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const double & q_i = q[iDimension];
    x_[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }
}


double ClassicSVfitIntegrand::EvalMET_TF(unsigned int iComponent) const{

  return EvalMET_TF(measuredMETx_[iComponent], measuredMETy_[iComponent],
                    covMET_[iComponent]);
}


double ClassicSVfitIntegrand::EvalMET_TF(const double & aMETx, const double & aMETy, const TMatrixD& covMET) const{

// determine transfer matrix for MET
    double invCovMETxx = covMET(1,1);
    double invCovMETxy = -covMET(0,1);
    double invCovMETyx = -covMET(1,0);
    double invCovMETyy = covMET(0,0);
    double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet)<1E-10){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
             <<"METx: "<<aMETy<<" METy: "<<aMETy
             << std::endl;
             //errorCode_ |= MatrixInversion; //FIXME violates const
    return 0;
  }
    double const_MET = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));

// evaluate transfer function for MET/hadronic recoil
  double residualX = aMETx - (nu1P4_.X() + nu2P4_.X());
  double residualY = aMETy - (nu1P4_.Y() + nu2P4_.Y());
#ifdef USE_SVFITTF
  if ( rhoHadTau_ != 0. ) {

    int tmpIndex = legIntegrationParams_[0].idx_VisPtShift_;
    double visPtShift1 = ( tmpIndex != -1 && !leg1isLep_ ) ? (1./x_[tmpIndex]) : 1.;

    tmpIndex = legIntegrationParams_[1].idx_VisPtShift_;
    double visPtShift2 = ( tmpIndex != -1 && !leg2isLep_ ) ? (1./x_[tmpIndex]) : 1.;
    if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 ) return 0.;

    residualX += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.px() + (visPtShift2 - 1.)*measuredTauLepton2_.px()));
    residualY += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.py() + (visPtShift2 - 1.)*measuredTauLepton2_.py()));
  }
#endif
  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
                 residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2/=covDet;

  if ( verbosity_ >= 2 ) {
      double sumNuPx = nu1P4_.X() + nu2P4_.X();
      double sumNuPy = nu1P4_.Y() + nu2P4_.Y();
    std::cout << "TF(met): recPx = " << aMETx << ", recPy = " << aMETy
                     << ", genPx = " << sumNuPx       << ", genPy = " << sumNuPy
                     << " pull2 = "<<pull2
                     << " prob = "<<const_MET*TMath::Exp(-0.5*pull2)
                     <<std::endl;
  }
  return const_MET*TMath::Exp(-0.5*pull2);
}

double
ClassicSVfitIntegrand::EvalPS(const double* q) const
{
  rescaleX(q);

  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand::Eval(const double*)>:" << std::endl;
    std::cout << " x = { ";
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << x_[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ != 0 ) { return 0.; }

  int tmpIndex = legIntegrationParams_[0].idx_VisPtShift_;
  double visPtShift1 = ( tmpIndex != -1 && !leg1isLep_ ) ? (1./x_[tmpIndex]) : 1.;

  tmpIndex = legIntegrationParams_[1].idx_VisPtShift_;
  double visPtShift2 = ( tmpIndex != -1 && !leg2isLep_ ) ? (1./x_[tmpIndex]) : 1.;
  if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 ) return 0.;

  //FIXME violates const computeVisMom(visPtShift1, visPtShift2);

  // compute visible energy fractions for both taus
  tmpIndex = legIntegrationParams_[0].idx_X_;
  assert(tmpIndex != -1);
  double x1_dash = x_[tmpIndex];
  double x1 = x1_dash/visPtShift1;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  double x2_dash = 0.0;
  tmpIndex = legIntegrationParams_[1].idx_X_;
  if (tmpIndex != -1) {
    x2_dash = x_[tmpIndex];
  }
  else {
    x2_dash = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).M2()/(diTauMassConstraint_ * diTauMassConstraint_)/x1_dash;
  }
  double x2 = x2_dash/visPtShift2;
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double nuEn = vis1En_*(1. - x1)/x1;

  tmpIndex = legIntegrationParams_[0].idx_mNuNu_;
  double nu1Mass = 0;
  double nu1P = nuEn;
  if(tmpIndex != -1){
    nu1Mass = TMath::Sqrt(x_[tmpIndex]);
    nu1P = TMath::Sqrt(TMath::Max(0., nuEn*nuEn - x_[tmpIndex]));
  }

  tmpIndex = legIntegrationParams_[0].idx_phi_;
  assert(tmpIndex != -1);
  double phiNu = x_[tmpIndex];
  double cosThetaNu = compCosThetaNuNu(vis1En_, vis1P_, leg1Mass2_, nuEn, nu1P, square(nu1Mass));
  if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) return 0.;

  double cosPhiNu, sinPhiNu, sinThetaNu;
  sincos(phiNu, &sinPhiNu, &cosPhiNu);
  double thetaNu = TMath::ACos(cosThetaNu);
  sinThetaNu = TMath::Sin(thetaNu);

  double nuPx_local = nu1P*cosPhiNu*sinThetaNu;
  double nuPy_local = nu1P*sinPhiNu*sinThetaNu;
  double nuPz_local = nu1P*cosThetaNu;
  double nuPx = nuPx_local*leg1eX_x_ + nuPy_local*leg1eY_x_ + nuPz_local*leg1eZ_x_;
  double nuPy = nuPx_local*leg1eX_y_ + nuPy_local*leg1eY_y_ + nuPz_local*leg1eZ_y_;
  double nuPz = nuPx_local*leg1eX_z_ + nuPy_local*leg1eY_z_ + nuPz_local*leg1eZ_z_;
  //std::cout << "nu1: En = " << nuEn << ", Pt = " << TMath::Sqrt(nu1Px*nu1Px + nu1Py*nu1Py) << std::endl;
  nu1P4_.SetPxPyPzE(nuPx, nuPy, nuPz, nuEn);

  double tau1En = vis1En_ + nuEn;
  double tau1Px = vis1P4_.px() + nuPx;
  double tau1Py = vis1P4_.py() + nuPy;
  double tau1Pz = vis1P4_.pz() + nuPz;
  //std::cout << "tau1: En = " << tau1En << ", Pt = " << TMath::Sqrt(tau1Px*tau1Px + tau1Py*tau1Py) << std::endl;
  tau1P4_.SetPxPyPzE(tau1Px, tau1Py, tau1Pz, tau1En);

  // compute neutrino and tau lepton four-vector for second tau
  nuEn = vis2En_*(1. - x2)/x2;

  tmpIndex = legIntegrationParams_[1].idx_mNuNu_;
  double nu2Mass = 0;
  double nu2P = nuEn;
  if(tmpIndex != -1){
    nu2Mass = TMath::Sqrt(x_[tmpIndex]);
    nu2P = TMath::Sqrt(TMath::Max(0., nuEn*nuEn - x_[tmpIndex]));
  }

  tmpIndex = legIntegrationParams_[1].idx_phi_;
  assert(tmpIndex != -1);
  phiNu = x_[tmpIndex];
  cosThetaNu = compCosThetaNuNu(vis2En_, vis2P_, leg2Mass2_, nuEn, nu2P, square(nu2Mass));
  if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) return 0.;

  sincos(phiNu, &sinPhiNu, &cosPhiNu);
  thetaNu = TMath::ACos(cosThetaNu);
  sinThetaNu = TMath::Sin(thetaNu);

  nuPx_local = nu2P*cosPhiNu*sinThetaNu;
  nuPy_local = nu2P*sinPhiNu*sinThetaNu;
  nuPz_local = nu2P*cosThetaNu;
  nuPx = nuPx_local*leg2eX_x_ + nuPy_local*leg2eY_x_ + nuPz_local*leg2eZ_x_;
  nuPy = nuPx_local*leg2eX_y_ + nuPy_local*leg2eY_y_ + nuPz_local*leg2eZ_y_;
  nuPz = nuPx_local*leg2eX_z_ + nuPy_local*leg2eY_z_ + nuPz_local*leg2eZ_z_;
  //std::cout << "nu2: En = " << nuEn << ", Pt = " << TMath::Sqrt(nuPx*nuPx + nuPy*nuPy) << std::endl;
  nu2P4_.SetPxPyPzE(nuPx, nuPy, nuPz, nuEn);

  double tau2En = vis2En_ + nuEn;
  double tau2Px = vis2P4_.px() + nuPx;
  double tau2Py = vis2P4_.py() + nuPy;
  double tau2Pz = vis2P4_.pz() + nuPz;
  //std::cout << "tau2: En = " << tau2En << ", Pt = " << TMath::Sqrt(tau2Px*tau2Px + tau2Py*tau2Py) << std::endl;
  tau2P4_.SetPxPyPzE(tau2Px, tau2Py, tau2Pz, tau2En);

  if ( verbosity_ >= 2 ) {
    std::cout << "leg1: En = " << vis1P4_.energy() << ", Px = " << vis1P4_.px() << ", Py = " << vis1P4_.py() << ", Pz = " << vis1P4_.pz() << ";"
              << " Pt = " << vis1P4_.pt() << ", eta = " << vis1P4_.eta() << ", phi = " << vis1P4_.phi() << ", mass = " << vis1P4_.mass()
              << " (x = " << x1 << ")" << std::endl;
    std::cout << "tau1: En = " << tau1P4_.energy() << ", Px = " << tau1P4_.px() << ", Py = " << tau1P4_.py() << ", Pz = " << tau1P4_.pz() << ";"
              << " Pt = " << tau1P4_.pt() << ", eta = " << tau1P4_.eta() << ", phi = " << tau1P4_.phi() << std::endl;
    std::cout << "nu1: En = " << nu1P4_.energy() << ", Px = " << nu1P4_.px() << ", Py = " << nu1P4_.py() << ", Pz = " << nu1P4_.pz() << ";"
              << " Pt = " << nu1P4_.pt() << ", eta = " << nu1P4_.eta() << ", phi = " << nu1P4_.phi() << ", mass = " << nu1P4_.mass() << std::endl;
    //double angle1 = compAngle(vis1P4_, nu1P4);
    //std::cout << "angle(vis1, nu1) = " << angle1 << std::endl;
    //double phiInvis1 = compPhiInvis(vis1P4_, nu1P4);
    //std::cout << "phiInvis1 = " << phiInvis1 << std::endl;
    std::cout << "leg2: En = " << vis2P4_.energy() << ", Px = " << vis2P4_.px() << ", Py = " << vis2P4_.py() << ", Pz = " << vis2P4_.pz() << ";"
              << " Pt = " << vis2P4_.pt() << ", eta = " << vis2P4_.eta() << ", phi = " << vis2P4_.phi() << ", mass = " << vis2P4_.mass()
              << " (x = " << x2 << ")" << std::endl;
    std::cout << "tau2: En = " << tau2P4_.energy() << ", Px = " << tau2P4_.px() << ", Py = " << tau2P4_.py() << ", Pz = " << tau2P4_.pz() << ";"
              << " Pt = " << tau2P4_.pt() << ", eta = " << tau2P4_.eta() << ", phi = " << tau2P4_.phi() << std::endl;
    std::cout << "nu2: En = " << nu2P4_.energy() << ", Px = " << nu2P4_.px() << ", Py = " << nu2P4_.py() << ", Pz = " << nu2P4_.pz() << ";"
              << " Pt = " << nu2P4_.pt() << ", eta = " << nu2P4_.eta() << ", phi = " << nu2P4_.phi() << ", mass = " << nu2P4_.mass() << std::endl;
    //double angle2 = compAngle(vis2P4_, nu2P4);
    //std::cout << "angle(vis2, nu2) = " << angle2 << std::endl;
    //double phiInvis2 = compPhiInvis(vis2P4_, nu2P4);
    //std::cout << "phiInvis2 = " << phiInvis2 << std::endl;
  }
  double prob_TF = 1.0;

#ifdef USE_SVFITTF
  // evaluate transfer functions for tau energy reconstruction
  if ( useHadTauTF_ && legIntegrationParams_[0].idx_VisPtShift_ != -1 && !leg1isLep_ ) {
    double prob_TF_leg1 = (*hadTauTF1_)(measuredTauLepton1_.pt(), vis1P4_.pt(), vis1P4_.eta());
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg1): recPt = " << measuredTauLepton1_.pt() << ", genPt = " << vis1P4_.pt()
                << ", genEta = " << vis1P4_.eta() << " --> prob = " << prob_TF_leg1 << std::endl;
    }
    prob_TF *= prob_TF_leg1;
  }
  if ( useHadTauTF_ && legIntegrationParams_[1].idx_VisPtShift_ != -1 && !leg2isLep_ ) {
    double prob_TF_leg2 = (*hadTauTF2_)(measuredTauLepton2_.pt(), vis2P4_.pt(), vis2P4_.eta());
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg2): recPt = " << measuredTauLepton2_.pt() << ", genPt = " << vis2P4_.pt()
                << ", genEta = " << vis2P4_.eta() << " --> prob = " << prob_TF_leg2 << std::endl;
    }
    prob_TF *= prob_TF_leg2;
  }
#endif

  double prob_PS_and_tauDecay = classic_svFit::constFactor;
  double prob_tauDecay_leg1 = 0.;
  if ( leg1isLep_ ) {
    prob_tauDecay_leg1 = compPSfactor_tauToLepDecay(x1, vis1P4_.E(), vis1P4_.P(), leg1Mass_, nu1P4_.E(), nu1P4_.P(), nu1Mass);
  } else {
    prob_tauDecay_leg1 = compPSfactor_tauToHadDecay(x1, vis1P4_.E(), vis1P4_.P(), leg1Mass_, nu1P4_.E(), nu1P4_.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg1;
  double prob_tauDecay_leg2 = 0.;
  if ( leg2isLep_ ) {
    prob_tauDecay_leg2 = compPSfactor_tauToLepDecay(x2, vis2P4_.E(), vis2P4_.P(), leg2Mass_, nu2P4_.E(), nu2P4_.P(), nu2Mass);
  } else {
    prob_tauDecay_leg2 = compPSfactor_tauToHadDecay(x2, vis2P4_.E(), vis2P4_.P(), leg2Mass_, nu2P4_.E(), nu2P4_.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg2;
  prob_PS_and_tauDecay *= classic_svFit::matrixElementNorm;

  double mTauTau = (tau1P4_ + tau2P4_).mass();
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
    jacobiFactor *= (2.0*x2/diTauMassConstraint_);
  }

  double prob = prob_PS_and_tauDecay*prob_TF*prob_logM*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "mTauTau = "<<mTauTau
              << " prob: PS+decay = " << prob_PS_and_tauDecay << ","
              << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
}
  if (TMath::IsNaN(prob) ) {
    //std::cerr << "Warning: prob = " << prob << " (PS+decay = " << prob_PS_and_tauDecay << ","
      //        << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor << ") --> setting prob = 0 !!" << std::endl;
    prob = 0.;
  }

  return prob;
}

double ClassicSVfitIntegrand::Eval(const double* x, unsigned int iComponent) const
{
  // CV: fill histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
    if(histogramAdapter_){
      histogramAdapter_->setTau1P4(tau1P4_);
      histogramAdapter_->setTau2P4(tau2P4_);
    }

    if(iComponent==0) phaseSpaceComponentCache_ = EvalPS(x);
    double metTF = EvalMET_TF(iComponent);

    return phaseSpaceComponentCache_*metTF;
}
