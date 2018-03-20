#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrandBase.h"

#include <TMath.h>
#include <TString.h> // Form
#include <Math/VectorUtil.h>

#include <math.h>

using namespace classic_svFit;

ClassicSVfitIntegrandBase::ClassicSVfitIntegrandBase(int verbosity)
  : numTaus_(0)
#ifdef USE_SVFITTF
  , useHadTauTF_(false)
  , rhoHadTau_(0.)
#endif
  , numDimensions_(0)
  , maxNumberOfDimensions_(0)
  , xMin_(nullptr)
  , xMax_(nullptr)
  , x_(nullptr)
  , addLogM_fixed_(false)
  , addLogM_fixed_power_(0.)
  , addLogM_dynamic_(false)
  , addLogM_dynamic_formula_(0)
  , errorCode_(0)
  , verbosity_(verbosity)
{}

ClassicSVfitIntegrandBase::~ClassicSVfitIntegrandBase()
{
#ifdef USE_SVFITTF
  for ( std::vector<const HadTauTFBase*>::iterator hadTauTF = hadTauTFs_.begin();
	hadTauTF != hadTauTFs_.end(); ++hadTauTF ) {
    delete (*hadTauTF);
  }
#endif

  delete xMin_;
  delete xMax_;
  delete x_;

  delete addLogM_dynamic_formula_;
}


void ClassicSVfitIntegrandBase::setVerbosity(int aVerbosity)
{
  verbosity_ = aVerbosity;
}

void ClassicSVfitIntegrandBase::addLogM_fixed(bool value, double power)
{
  addLogM_fixed_ = value;
  addLogM_fixed_power_ = power;
  if ( addLogM_fixed_ && addLogM_dynamic_ ) {
    std::cerr << "Warning: simultaneous use of fixed and dynamic logM terms not supported --> disabling dynamic logM term !!" << std::endl;
    addLogM_dynamic_ = false;
  }
}

void ClassicSVfitIntegrandBase::addLogM_dynamic(bool value, const std::string& power)
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

void ClassicSVfitIntegrandBase::setLegIntegrationParams(unsigned int iLeg, const classic_svFit::integrationParameters& aParams)
{ 
  assert(iLeg < legIntegrationParams_.size());
  legIntegrationParams_[iLeg] = aParams;
}

void ClassicSVfitIntegrandBase::setNumDimensions(unsigned numDimensions) 
{ 
  assert(numDimensions <= maxNumberOfDimensions_);
  numDimensions_ = numDimensions; 
}

void ClassicSVfitIntegrandBase::setIntegrationRanges(const double* xl, const double* xh)
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xh[iDimension];
  }
}

#ifdef USE_SVFITTF
void ClassicSVfitIntegrandBase::setHadTauTF(const HadTauTFBase* hadTauTF)
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

void ClassicSVfitIntegrandBase::enableHadTauTF()
{
  if ( !(hadTauTFs_.size() == numTaus_) ) {
    std::cerr << "No tau pT transfer functions defined, call 'setHadTauTF' function first !!" << std::endl;
    assert(0);
  }
  useHadTauTF_ = true;
}

void ClassicSVfitIntegrandBase::disableHadTauTF()
{
  useHadTauTF_ = false;
}

void ClassicSVfitIntegrandBase::setRhoHadTau(double rhoHadTau)
{
  rhoHadTau_ = rhoHadTau;
}
#endif


void ClassicSVfitIntegrandBase::setLeptonInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  // reset 'LeptonNumber' and 'MatrixInversion' error codes
  errorCode_ &= (errorCode_ ^ LeptonNumber);
  errorCode_ &= (errorCode_ ^ MatrixInversion);

  if ( measuredTauLeptons.size() != numTaus_ ) {
    std::cerr << "Error: Number of MeasuredTauLeptons is not equal to " << numTaus_ << " !!" << std::endl;
    errorCode_ |= LeptonNumber;
  }

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

void ClassicSVfitIntegrandBase::addMETEstimate(double measuredMETx, double measuredMETy, const TMatrixD& covMET)
{
  measuredMETx_.push_back(measuredMETx);
  measuredMETy_.push_back(measuredMETy);
  covMET_.push_back(covMET);
}

int ClassicSVfitIntegrandBase::getMETComponentsSize() const 
{
  return measuredMETx_.size();
}

void ClassicSVfitIntegrandBase::clearMET()
{
  measuredMETx_.clear();
  measuredMETy_.clear();
  covMET_.clear();
}

void ClassicSVfitIntegrandBase::rescaleX(const double* q) const
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const double& q_i = q[iDimension];
    x_[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }
}

double ClassicSVfitIntegrandBase::EvalMET_TF(unsigned int iComponent) const
{
  return EvalMET_TF(measuredMETx_[iComponent], measuredMETy_[iComponent], covMET_[iComponent]);
}

double ClassicSVfitIntegrandBase::EvalMET_TF(double aMETx, double aMETy, const TMatrixD& covMET) const
{
  // determine transfer matrix for MET
  double invCovMETxx =  covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy =  covMET(0,0);
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


