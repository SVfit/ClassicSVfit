#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"

#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include <TMath.h>
#include <TString.h> // Form
#include <Math/VectorUtil.h>

#include <math.h>

using namespace classic_svFit;

/// global function pointer, needed for Markov Chain integration
const ClassicSVfitIntegrand* ClassicSVfitIntegrand::gSVfitIntegrand = 0;

ClassicSVfitIntegrand::ClassicSVfitIntegrand(int verbosity)
  : numTaus_(2)
  , fittedTauLepton1_(0, verbosity)
  , fittedTauLepton2_(1, verbosity)
#ifdef USE_SVFITTF
  , useHadTauTF_(false)
  , rhoHadTau_(0.)
#endif
  , numDimensions_(0)
  , addLogM_fixed_(true)
  , addLogM_fixed_power_(6.) // CV: best compatibility with "old" SVfitStandalone algorithm
  , addLogM_dynamic_(false)
  , addLogM_dynamic_formula_(0)
  , diTauMassConstraint_(-1.)
  , errorCode_(0)
  , histogramAdapter_(0)
  , verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<ClassicSVfitIntegrand::ClassicSVfitIntegrand>:" << std::endl;
  }

  fittedTauLeptons_.resize(numTaus_);
  fittedTauLeptons_[0] = &fittedTauLepton1_;
  fittedTauLeptons_[1] = &fittedTauLepton2_;

  // set global function pointer to this
  gSVfitIntegrand = this;
}

ClassicSVfitIntegrand::~ClassicSVfitIntegrand()
{
  if ( verbosity_ ) {
    std::cout << "<ClassicSVfitIntegrand::~ClassicSVfitIntegrand>:" << std::endl;
  }

#ifdef USE_SVFITTF
  for ( std::vector<const HadTauTFBase*>::iterator hadTauTF = hadTauTFs_.begin();
	hadTauTF != hadTauTFs_.end(); ++hadTauTF ) {
    delete (*hadTauTF);
  }
#endif

  delete addLogM_dynamic_formula_;
}


void ClassicSVfitIntegrand::setVerbosity(int aVerbosity)
{
  verbosity_ = aVerbosity;
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
  diTauMassConstraint2_ = square(diTauMassConstraint_);
}

void ClassicSVfitIntegrand::setHistogramAdapter(HistogramAdapterDiTau* histogramAdapter)
{
  histogramAdapter_ = histogramAdapter;
}

void ClassicSVfitIntegrand::setLegIntegrationParams(unsigned int iLeg, const classic_svFit::integrationParameters& aParams)
{ 
  legIntegrationParams_[iLeg] = aParams;
}

void ClassicSVfitIntegrand::setNumDimensions(unsigned numDimensions) 
{ 
  numDimensions_ = numDimensions; 
}

void ClassicSVfitIntegrand::setIntegrationRanges(const double* xl, const double* xu)
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
  }
}

#ifdef USE_SVFITTF
void ClassicSVfitIntegrand::setHadTauTF(const HadTauTFBase* hadTauTF)
{
  for ( std::vector<const HadTauTFBase*>::iterator hadTauTF = hadTauTFs_.begin();
	hadTauTF != hadTauTFs_.end(); ++hadTauTF ) {
    delete (*hadTauTF);
  }
  hadTauTFs_.clear();
  for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
    hadTauTFs_.push_back(hadTauTF->Clone(Form("leg%i" + iTau)));
  }
}

void ClassicSVfitIntegrand::enableHadTauTF()
{
  if ( !(hadTauTFs_.size() == numTaus_) ) {
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


void ClassicSVfitIntegrand::setLeptonInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand::setLeptonInputs>:" << std::endl;
  }

  // reset 'LeptonNumber' and 'MatrixInversion' error codes
  errorCode_ &= (errorCode_ ^ LeptonNumber);
  errorCode_ &= (errorCode_ ^ MatrixInversion);

  if ( measuredTauLeptons.size() != numTaus_ ) {
    std::cerr << "Error: Number of MeasuredTauLeptons is not equal to " << numTaus_ << " !!" << std::endl;
    errorCode_ |= LeptonNumber;
  }
  
  // set momenta of visible tau decay products, reset momenta of reconstructed tau leptons
  measuredTauLepton1_ = measuredTauLeptons[0];
  fittedTauLepton1_.setMeasuredTauLepton(measuredTauLepton1_);
  leg1isLeptonicTauDecay_ = measuredTauLepton1_.isLeptonicTauDecay();
  leg1isHadronicTauDecay_ = measuredTauLepton1_.isHadronicTauDecay();
  leg1isPrompt_ = measuredTauLepton1_.isPrompt();
  measuredTauLepton2_ = measuredTauLeptons[1];
  fittedTauLepton2_.setMeasuredTauLepton(measuredTauLepton2_);
  leg2isLeptonicTauDecay_ = measuredTauLepton2_.isLeptonicTauDecay();
  leg2isHadronicTauDecay_ = measuredTauLepton2_.isHadronicTauDecay();
  leg2isPrompt_ = measuredTauLepton2_.isPrompt();

  mVis_measured_ = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).mass();
  if ( verbosity_ >= 2 ) {
    std::cout << "mVis = " << mVis_measured_ << std::endl;
  }
  mVis2_measured_ = square(mVis_measured_);

  phaseSpaceComponentCache_ = 0;

#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) {
    for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
      const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
      const MeasuredTauLepton& measuredTauLepton = fittedTauLepton->getMeasuredTauLepton();
      if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay ) {
	hadTauTFs_[iTau]->setDecayMode(measuredTauLepton.decayMode());
      }
    }
  }
#endif
}

void ClassicSVfitIntegrand::addMETEstimate(double measuredMETx, double measuredMETy, const TMatrixD& covMET)
{
  measuredMETx_.push_back(measuredMETx);
  measuredMETy_.push_back(measuredMETy);
  covMET_.push_back(covMET);
}

int ClassicSVfitIntegrand::getMETComponentsSize() const 
{
  return measuredMETx_.size();
}

void ClassicSVfitIntegrand::clearMET()
{
  measuredMETx_.clear();
  measuredMETy_.clear();
  covMET_.clear();
}

void ClassicSVfitIntegrand::rescaleX(const double* q) const
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const double& q_i = q[iDimension];
    x_[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }
}

double ClassicSVfitIntegrand::EvalMET_TF(unsigned int iComponent) const
{
  return EvalMET_TF(measuredMETx_[iComponent], measuredMETy_[iComponent], covMET_[iComponent]);
}

double ClassicSVfitIntegrand::EvalMET_TF(double aMETx, double aMETy, const TMatrixD& covMET) const
{
  // determine transfer matrix for MET
  double invCovMETxx = covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy = covMET(0,0);
  double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet) < 1.e-10 ){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!" << std::endl;
    errorCode_ |= MatrixInversion;
    return 0;
  }
  double const_MET = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));

  // compute sum of momenta of all neutrinos produced in tau decays
  double sumNuPx = 0.;
  double sumNuPy = 0.;
  for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
    const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
    sumNuPx += fittedTauLepton->nuP4().px();
    sumNuPy += fittedTauLepton->nuP4().py();
  }

  // evaluate transfer function for MET/hadronic recoil
  double residualX = aMETx - sumNuPx;
  double residualY = aMETy - sumNuPy;
#ifdef USE_SVFITTF
  if ( rhoHadTau_ != 0. ) {
    for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
      const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
      const MeasuredTauLepton& measuredTauLepton = fittedTauLepton->getMeasuredTauLepton();
      if ( measuredTauLepton.isHadronicTauDecay() ) {
	int idx_visPtShift = legIntegrationParams_[iTau].idx_VisPtShift_;
	if ( idx_visPtShift != -1 ) {
	  double visPtShift = 1./x_[idx_visPtShift];
	  if ( visPtShift < 1.e-2 ) continue;
	  residualX += (rhoHadTau_*(visPtShift - 1.)*measuredTauLepton.px());
	  residualY += (rhoHadTau_*(visPtShift - 1.)*measuredTauLepton.py());
	}
      }
    }
  }
#endif
  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
                 residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2 /= covDet;
  double prob = const_MET*TMath::Exp(-0.5*pull2);

  if ( verbosity_ >= 2 ) {    
    std::cout << "TF(met): recPx = " << aMETx << ", recPy = " << aMETy << ","
	      << " genPx = " << sumNuPx << ", genPy = " << sumNuPy << ","
	      << " pull2 = " << pull2 << ", prob = " << prob << std::endl;
  }
  return prob;
}

double
ClassicSVfitIntegrand::EvalPS(const double* q) const
{
  rescaleX(q);

  if ( verbosity_ >= 2 ) {
    std::cout << "<ClassicSVfitIntegrand::EvalPS(const double*)>:" << std::endl;
    std::cout << " x = { ";
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << x_[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ & MatrixInversion ||
       errorCode_ & LeptonNumber    ||
       errorCode_ & TestMass        ) {
    return 0.; 
  }

  int idx_visPtShift1 = legIntegrationParams_[0].idx_VisPtShift_;
  double visPtShift1 = ( useHadTauTF_ && idx_visPtShift1 != -1 && !leg1isLeptonicTauDecay_ ) ? (1./x_[idx_visPtShift1]) : 1.;
  int idx_visPtShift2 = legIntegrationParams_[1].idx_VisPtShift_;
  double visPtShift2 = ( useHadTauTF_ && idx_visPtShift2 != -1 && !leg2isLeptonicTauDecay_ ) ? (1./x_[idx_visPtShift2]) : 1.;
  if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 ) return 0.;

  // scale momenta of visible tau decays products
  fittedTauLepton1_.updateVisMomentum(visPtShift1);
  fittedTauLepton2_.updateVisMomentum(visPtShift2);

  // compute visible energy fractions for both taus
  double x1_dash = 1.;
  if ( !leg1isPrompt_ ) {
    int idx_x1 = legIntegrationParams_[0].idx_X_;
    assert(idx_x1 != -1);
    x1_dash = x_[idx_x1];
  }
  double x1 = x1_dash/visPtShift1;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  double x2_dash = 1.;
  if ( !leg2isPrompt_ ) {
    int idx_x2 = legIntegrationParams_[1].idx_X_;
    if ( idx_x2 != -1 ) {
      x2_dash = x_[idx_x2];
    } else {
      x2_dash = (mVis2_measured_/diTauMassConstraint2_)/x1_dash;
    }
  }
  double x2 = x2_dash/visPtShift2;
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  // compute neutrino and tau lepton momenta 
  if ( !leg1isPrompt_ ) {
    int idx_phiNu1 = legIntegrationParams_[0].idx_phi_;
    assert(idx_phiNu1 != -1);
    double phiNu1 = x_[idx_phiNu1];
    int idx_nu1Mass = legIntegrationParams_[0].idx_mNuNu_;
    double nu1Mass = ( idx_nu1Mass != -1 ) ? TMath::Sqrt(x_[idx_nu1Mass]) : 0.;
    fittedTauLepton1_.updateTauMomentum(x1, phiNu1, nu1Mass);
    //std::cout << "fittedTauLepton1: errorCode = " << fittedTauLepton1_.errorCode() << std::endl;
    if ( fittedTauLepton1_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( !leg2isPrompt_ ) {
    int idx_phiNu2 = legIntegrationParams_[1].idx_phi_;
    assert(idx_phiNu2 != -1);
    double phiNu2 = x_[idx_phiNu2];
    int idx_nu2Mass = legIntegrationParams_[1].idx_mNuNu_;
    double nu2Mass = ( idx_nu2Mass != -1 ) ? TMath::Sqrt(x_[idx_nu2Mass]) : 0.;
    fittedTauLepton2_.updateTauMomentum(x2, phiNu2, nu2Mass);
    //std::cout << "fittedTauLepton2: errorCode = " << fittedTauLepton1_.errorCode() << std::endl;
    if ( fittedTauLepton2_.errorCode() != FittedTauLepton::None ) {
      errorCode_ |= TauDecayParameters;
      return 0.;
    }
  }

  if ( verbosity_ >= 2 ) {
    for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
      const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
      const LorentzVector& visP4 = fittedTauLepton->visP4();
      const LorentzVector& nuP4 = fittedTauLepton->nuP4();
      const LorentzVector& tauP4 = fittedTauLepton->tauP4();
      std::cout << "leg" << (iTau + 1) << ": En = " << visP4.E() << ", Px = " << visP4.px()
		<< ", Py = " << visP4.py() << ", Pz = " << visP4.pz() << ";"
		<< " Pt = " << visP4.pt() << ", eta = " << visP4.eta()
		<< ", phi = " << visP4.phi() << ", mass = " << visP4.mass()
		<< " (x = " << fittedTauLepton->x() << ")" << std::endl;
      std::cout << "tau" << (iTau + 1) << ": En = " << tauP4.E() << ", Px = " << tauP4.px() << ", Py = " << tauP4.py() << ", Pz = " << tauP4.pz() << ";"
		<< " Pt = " << tauP4.pt() << ", eta = " << tauP4.eta() << ", phi = " << tauP4.phi() << std::endl;
      std::cout << "nu" << (iTau + 1) << ": En = " << nuP4.E() << ", Px = " << nuP4.px() << ", Py = " << nuP4.py() << ", Pz = " << nuP4.pz() << ";"
		<< " Pt = " << nuP4.pt() << ", eta = " << nuP4.eta() << ", phi = " << nuP4.phi() << ", mass = " << nuP4.mass() << std::endl;
    }
  }

  double prob_PS_and_tauDecay = classic_svFit::constFactor;
  double prob_tauDecay = 1.;
  double prob_TF = 1.;
  for ( unsigned iTau = 0; iTau < numTaus_; ++iTau ) {
    const FittedTauLepton* fittedTauLepton = fittedTauLeptons_[iTau];
    const MeasuredTauLepton& measuredTauLepton = fittedTauLepton->getMeasuredTauLepton();
    double x = fittedTauLepton->x();
    double nuMass = fittedTauLepton->nuMass();
    const LorentzVector& visP4 = fittedTauLepton->visP4();
    const LorentzVector& nuP4 = fittedTauLepton->nuP4();

    // evaluate tau decay matrix elements
    double prob = 1.;
    if      ( measuredTauLepton.isLeptonicTauDecay() ) prob = compPSfactor_tauToLepDecay(x, visP4.E(), visP4.P(), measuredTauLepton.mass(), nuP4.E(), nuP4.P(), nuMass);
    else if ( measuredTauLepton.isHadronicTauDecay() ) prob = compPSfactor_tauToHadDecay(x, visP4.E(), visP4.P(), measuredTauLepton.mass(), nuP4.E(), nuP4.P());
    prob_tauDecay *= prob;

    // evaluate transfer functions for tau energy reconstruction
#ifdef USE_SVFITTF
    if ( useHadTauTF_ && legIntegrationParams_[iTau].idx_VisPtShift_ != -1 && measuredTauLepton.isHadronicTauDecay() ) {
      double prob = (*hadTauTFs_[iTau])(measuredTauLepton.pt(), visP4.pt(), visP4.eta());
      if ( verbosity_ >= 2 ) {
	std::cout << "TF(leg" << iTau << "): recPt = " << measuredTauLepton.pt() << ", genPt = " << visP4.pt()
		  << ", genEta = " << visP4.eta() << " --> prob = " << prob << std::endl;
      }
      prob_TF *= prob;
    }
#endif
  }
  prob_PS_and_tauDecay *= prob_tauDecay;
  prob_PS_and_tauDecay *= classic_svFit::matrixElementNorm;

  double mTauTau = (fittedTauLepton1_.tauP4() + fittedTauLepton2_.tauP4()).mass();
  double prob_logM = 1.;
  if ( addLogM_fixed_ ) {
    prob_logM = 1./TMath::Power(TMath::Max(1., mTauTau), addLogM_fixed_power_);
  }
  if ( addLogM_dynamic_ ) {
    double addLogM_power = addLogM_dynamic_formula_->Eval(mTauTau);
    prob_logM = 1./TMath::Power(TMath::Max(1., mTauTau), TMath::Max(0., addLogM_power));
  }

  double jacobiFactor = 1./(visPtShift1*visPtShift2); // product of derrivatives dx1/dx1' and dx2/dx2' for parametrization of x1, x2 by x1', x2'
  if ( diTauMassConstraint_ > 0. ) {
    jacobiFactor *= (2.*x2/diTauMassConstraint_);
  }

  //static int numCalls = 0;
  //++numCalls;
  //std::cout << "call #" << numCalls << ":" << std::endl;
  //if ( numCalls > 100 ) assert(0);

  double prob = prob_PS_and_tauDecay*prob_TF*prob_logM*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "mTauTau = " << mTauTau << std::endl;
    std::cout << "prob: PS+decay = " << prob_PS_and_tauDecay << ","
              << " TF = " << prob_TF << ", log(M) = " << prob_logM << ", Jacobi = " << jacobiFactor 
	      << " --> returning " << prob << std::endl;
  }
  if ( TMath::IsNaN(prob) ) {
    prob = 0.;
  }

  return prob;
}

double ClassicSVfitIntegrand::Eval(const double* x, unsigned int iComponent) const
{
  if ( iComponent == 0 ) {
    phaseSpaceComponentCache_ = EvalPS(x);
  }
  if ( phaseSpaceComponentCache_ < 1.e-300 ) return 0.;
  double prob = phaseSpaceComponentCache_*EvalMET_TF(iComponent);
  if ( verbosity_ >= 2 ) {
    std::cout<<" metTF: " << EvalMET_TF(iComponent)
	     << " phaseSpaceComponentCache_: " << phaseSpaceComponentCache_
	     << " --> returning " << prob << std::endl;
  }
  if ( histogramAdapter_ && prob > 1.e-300 ){
    histogramAdapter_->setTau1And2P4(fittedTauLepton1_.tauP4(), fittedTauLepton2_.tauP4());
  }
  return prob;
}
