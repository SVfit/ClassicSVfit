#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitBase.h"

#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include <TGraphErrors.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace classic_svFit;

ClassicSVfitBase::ClassicSVfitBase(int verbosity)
  : integrand_(0)
  , intAlgo_(0)
  , maxObjFunctionCalls_(100000)
  , treeFileName_("")
  , likelihoodFileName_("")
  , numDimensions_(0)
  , xl_(nullptr)
  , xh_(nullptr)
  , isValidSolution_(false)
  , useHadTauTF_(false)
  , clock_(nullptr)
  , numSeconds_cpu_(-1.)
  , numSeconds_real_(-1.)
  , verbosity_(verbosity)
{
  clock_ = new TBenchmark();
}

ClassicSVfitBase::~ClassicSVfitBase()
{
  delete integrand_;

  if ( intAlgo_ ) {
    delete intAlgo_;
  }

  delete [] xl_;
  delete [] xh_;

  delete clock_;
}

void ClassicSVfitBase::setVerbosity(int aVerbosity)
{
  verbosity_ = aVerbosity;
  integrand_->setVerbosity(verbosity_);
}

void ClassicSVfitBase::addLogM_fixed(bool value, double power)
{
  integrand_->addLogM_fixed(value, power);
}

void ClassicSVfitBase::addLogM_dynamic(bool value, const std::string& power)
{
  integrand_->addLogM_dynamic(value, power);
}

#ifdef USE_SVFITTF
void ClassicSVfitBase::setHadTauTF(const HadTauTFBase* hadTauTF)
{
  integrand_->setHadTauTF(hadTauTF);
}

void ClassicSVfitBase::enableHadTauTF()
{
  integrand_->enableHadTauTF();
  useHadTauTF_ = true;
}

void ClassicSVfitBase::disableHadTauTF()
{
  integrand_->disableHadTauTF();
  useHadTauTF_ = false;
}

void ClassicSVfitBase::setRhoHadTau(double rhoHadTau)
{
  integrand_->setRhoHadTau(rhoHadTau);
}
#endif


void ClassicSVfitBase::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  maxObjFunctionCalls_ = maxObjFunctionCalls;
}

void ClassicSVfitBase::setLikelihoodFileName(const std::string& likelihoodFileName)
{
  likelihoodFileName_ = likelihoodFileName;
}

void ClassicSVfitBase::setTreeFileName(const std::string& treeFileName)
{
  treeFileName_ = treeFileName;
}

bool ClassicSVfitBase::isValidSolution() const 
{
  return isValidSolution_;
}

double ClassicSVfitBase::getComputingTime_cpu() const 
{
  return numSeconds_cpu_;
}

double ClassicSVfitBase::getComputingTime_real() const 
{
  return numSeconds_real_;
}

void ClassicSVfitBase::initializeMCIntegrator()
{
  //unsigned numChains = TMath::Nint(maxObjFunctionCalls_/100000.);
  unsigned numChains = 1;
  unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls_/numChains);
  unsigned numIterSampling = TMath::Nint(0.90*maxObjFunctionCalls_/numChains);
  unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.20*numIterBurnin);
  unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.60*numIterBurnin);
  if ( treeFileName_ == "" && verbosity_ >= 2 ) {
    treeFileName_ = "SVfitIntegratorMarkovChain_ClassicSVfit.root";
  }
  intAlgo_ = new SVfitIntegratorMarkovChain(
    "uniform",
    numIterBurnin, numIterSampling, numIterSimAnnealingPhase1, numIterSimAnnealingPhase2,
    15., 1. - 1./(0.1*numIterBurnin),
    numChains, 100,
    1.e-2, 0.71,
    treeFileName_.data(),
    0);
}

void ClassicSVfitBase::printMET(double measuredMETx, double measuredMETy, const TMatrixD& covMET) const
{
  std::cout << "MET: Px = " << measuredMETx << ", Py = " <<  measuredMETy<< std::endl;
  std::cout << "covMET:" << std::endl;
  covMET.Print();
  TMatrixDSym covMET_sym(2);
  covMET_sym(0,0) = covMET[0][0];
  covMET_sym(0,1) = covMET[0][1];
  covMET_sym(1,0) = covMET[1][0];
  covMET_sym(1,1) = covMET[1][1];
  TMatrixD EigenVectors(2,2);
  EigenVectors = TMatrixDSymEigen(covMET_sym).GetEigenVectors();
  std::cout << "Eigenvectors =  { " << EigenVectors(0,0) << ", " << EigenVectors(1,0) << " (phi = " << TMath::ATan2(EigenVectors(1,0), EigenVectors(0,0)) << ") },"
	    << " { " << EigenVectors(0,1) << ", " << EigenVectors(1,1) << " (phi = " << TMath::ATan2(EigenVectors(1,1), EigenVectors(0,1)) << ") }" << std::endl;
  TVectorD EigenValues(2);
  EigenValues = TMatrixDSymEigen(covMET_sym).GetEigenValues();
  EigenValues(0) = TMath::Sqrt(EigenValues(0));
  EigenValues(1) = TMath::Sqrt(EigenValues(1));
  std::cout << "Eigenvalues = " << EigenValues(0) << ", " << EigenValues(1) << std::endl;
}

void ClassicSVfitBase::printLeptons() const 
{
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
    std::cout << "measuredTauLepton #" << idx << " (type = " << measuredTauLepton.type() << "): Pt = " << measuredTauLepton.pt() << ","
	      << " eta = " << measuredTauLepton.eta() << " (theta = " << measuredTauLepton.p3().theta() << ")" << ", phi = " << measuredTauLepton.phi() << ","
	      << " mass = " << measuredTauLepton.mass() << std::endl;
  }
}

void ClassicSVfitBase::printIntegrationRange() const 
{
  std::cout << "numDimensions = " << numDimensions_ << std::endl;

  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    std::cout << " fitParameter #" << iDimension << ": xl = " << xl_[iDimension] << ", xh = " << xh_[iDimension];
    for ( unsigned iLeg = 0; iLeg < legIntegrationParams_.size(); ++iLeg ) {
      if ( (int)iDimension == legIntegrationParams_[iLeg].idx_X_          ) std::cout << " (leg" << (iLeg + 1) << ":X)";
      if ( (int)iDimension == legIntegrationParams_[iLeg].idx_phi_        ) std::cout << " (leg" << (iLeg + 1) << ":phi)";
      if ( (int)iDimension == legIntegrationParams_[iLeg].idx_VisPtShift_ ) std::cout << " (leg" << (iLeg + 1) << ":VisPtShift)";
      if ( (int)iDimension == legIntegrationParams_[iLeg].idx_mNuNu_      ) std::cout << " (leg" << (iLeg + 1) << ":mNuNu)";
    }
    std::cout << std::endl;
  }
}

void ClassicSVfitBase::setLegIntegrationParams(unsigned iLeg, bool useMassConstraint)
{
  assert(iLeg < measuredTauLeptons_.size());
  const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[iLeg];
  
  if ( measuredTauLepton.type() == MeasuredTauLepton::kPrompt ) {
    if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay) {
      if ( useHadTauTF_ ) legIntegrationParams_[iLeg].idx_VisPtShift_ = numDimensions_++;
    }
  } else {
    if ( !useMassConstraint ) {
      legIntegrationParams_[iLeg].idx_X_ = numDimensions_++;
    }
    
    legIntegrationParams_[iLeg].idx_phi_ = numDimensions_++;
    
    if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay) {
      if ( useHadTauTF_ ) legIntegrationParams_[iLeg].idx_VisPtShift_ = numDimensions_++;
    } else {
      legIntegrationParams_[iLeg].idx_mNuNu_ = numDimensions_++;
    }
  }
  setIntegrationRanges(iLeg);
}

void ClassicSVfitBase::setIntegrationRanges(unsigned iLeg)
{
  const classic_svFit::integrationParameters& aIntParams = legIntegrationParams_[iLeg];

  if ( aIntParams.idx_X_ != -1 ) {
    xl_[aIntParams.idx_X_] = 0.;
#ifdef USE_SVFITTF
    if ( aIntParams.idx_VisPtShift_ != -1 ) {
      xh_[aIntParams.idx_X_] = 2.; // upper integration bound for x1' = visPtShift1*x1
    } else {
      xh_[aIntParams.idx_X_] = 1.;
    }
#else
    xh_[aIntParams.idx_X_] = 1.;
#endif
  }
  if ( aIntParams.idx_phi_ != -1 ) {
    xl_[aIntParams.idx_phi_] = -TMath::Pi();
    xh_[aIntParams.idx_phi_] = +TMath::Pi();
  }
  if ( aIntParams.idx_VisPtShift_ != -1 ) {
    xl_[aIntParams.idx_VisPtShift_] = 0.;
    xh_[aIntParams.idx_VisPtShift_] = 2.;
  }
  if ( aIntParams.idx_mNuNu_ != -1 ) {
    xl_[aIntParams.idx_mNuNu_] = 0.;
    xh_[aIntParams.idx_mNuNu_] = tauLeptonMass2;
  }
}

void ClassicSVfitBase::clearMET()
{ 
  integrand_->clearMET();
}

void ClassicSVfitBase::addMETEstimate(double measuredMETx, double measuredMETy, const TMatrixD& covMET)
{
  double metX = roundToNdigits(measuredMETx);
  double metY = roundToNdigits(measuredMETy);

  TMatrixD aCovMET(2,2);
  aCovMET[0][0] = roundToNdigits(covMET[0][0]);
  aCovMET[1][0] = roundToNdigits(covMET[1][0]);
  aCovMET[0][1] = roundToNdigits(covMET[0][1]);
  aCovMET[1][1] = roundToNdigits(covMET[1][1]);
  
  if ( verbosity_ >= 1 ) printMET(metX, metY, aCovMET);
  integrand_->addMETEstimate(metX, metY, aCovMET);
}
