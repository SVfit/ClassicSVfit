#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"

#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include <TGraphErrors.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>

using namespace classic_svFit;

namespace
{
  double g_C(double* x, size_t dim, void* param)
  {
    //std::cout << "<g_C>:" << std::endl;
    double retVal = ClassicSVfitIntegrand::gSVfitIntegrand->Eval(x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }
}

ClassicSVfit::ClassicSVfit(int verbosity)
  : integrand_(0),
    maxObjFunctionCalls_(100000),
    precision_(1.e-3),
    treeFileName_(""),
    numDimensions_(0),
    xl_(0),
    xu_(0),
    histogramAdapter_(new DiTauSystemHistogramAdapter()),
    likelihoodFileName_(""),
    isValidSolution_(false),
    useHadTauTF_(false),
    clock_(0),
    numSeconds_cpu_(-1.),
    numSeconds_real_(-1.),
    verbosity_(verbosity)
{
  integrand_ = new ClassicSVfitIntegrand(verbosity_);

  clock_ = new TBenchmark();
}

ClassicSVfit::~ClassicSVfit()
{
  delete histogramAdapter_;

  delete integrand_;

  delete clock_;
}

void ClassicSVfit::addLogM_fixed(bool value, double power)
{
  integrand_->addLogM_fixed(value, power);
}
void ClassicSVfit::addLogM_dynamic(bool value, const std::string& power)
{
  integrand_->addLogM_dynamic(value, power);
}

void ClassicSVfit::setDiTauMassConstraint(double diTauMass)
{
  diTauMassConstraint_ = diTauMass;
  integrand_->setDiTauMassConstraint(diTauMassConstraint_);
}

void ClassicSVfit::setTau1Constraint(bool tau1Constraint)
{
   tau1Constraint_ = tau1Constraint;
   integrand_->setTau1Constraint(tau1Constraint_);
}

#ifdef USE_SVFITTF
void ClassicSVfit::setHadTauTF(const HadTauTFBase* hadTauTF)
{
  integrand_->setHadTauTF(hadTauTF);
}
void ClassicSVfit::enableHadTauTF()
{
  integrand_->enableHadTauTF();
  useHadTauTF_ = true;
}
void ClassicSVfit::disableHadTauTF()
{
  integrand_->disableHadTauTF();
  useHadTauTF_ = false;
}

void ClassicSVfit::setRhoHadTau(double rhoHadTau)
{
  integrand_->setRhoHadTau(rhoHadTau);
}
#endif

void ClassicSVfit::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  maxObjFunctionCalls_ = maxObjFunctionCalls;
}

void ClassicSVfit::setLikelihoodFileName(const std::string& likelihoodFileName)
{
  likelihoodFileName_ = likelihoodFileName;
}

void ClassicSVfit::setTreeFileName(const std::string& treeFileName)
{
  treeFileName_ = treeFileName;
}

bool ClassicSVfit::isValidSolution() const { return isValidSolution_; }

double ClassicSVfit::getComputingTime_cpu() const { return numSeconds_cpu_; }
double ClassicSVfit::getComputingTime_real() const { return numSeconds_real_; }

namespace
{
  struct sortMeasuredTauLeptons
  {
    bool operator() (const MeasuredTauLepton& measuredTauLepton1, const MeasuredTauLepton& measuredTauLepton2)
    {
      if ( (measuredTauLepton1.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton1.type() == MeasuredTauLepton::kTauToMuDecay) &&
           measuredTauLepton2.type() == MeasuredTauLepton::kTauToHadDecay  ) return true;
      if ( (measuredTauLepton2.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2.type() == MeasuredTauLepton::kTauToMuDecay) &&
           measuredTauLepton1.type() == MeasuredTauLepton::kTauToHadDecay ) return false;
      return ( measuredTauLepton1.pt() > measuredTauLepton2.pt() );
    }
  };
}

void
ClassicSVfit::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<ClassicSVfit::integrate>:" << std::endl;
  }
  clock_->Reset();
  clock_->Start("<ClassicSVfit::integrate>");

  std::vector<MeasuredTauLepton> measuredTauLeptons_rounded;
  for ( std::vector<MeasuredTauLepton>::const_iterator measuredTauLepton = measuredTauLeptons.begin();
  measuredTauLepton != measuredTauLeptons.end(); ++measuredTauLepton ) {
    MeasuredTauLepton measuredTauLepton_rounded(
      measuredTauLepton->type(),
      roundToNdigits(measuredTauLepton->pt()),
      roundToNdigits(measuredTauLepton->eta()),
      roundToNdigits(measuredTauLepton->phi()),
      roundToNdigits(measuredTauLepton->mass()),
      measuredTauLepton->decayMode());
    measuredTauLeptons_rounded.push_back(measuredTauLepton_rounded);
  }
  std::sort(measuredTauLeptons_rounded.begin(), measuredTauLeptons_rounded.end(), sortMeasuredTauLeptons());
  measuredTauLeptons_ = measuredTauLeptons_rounded;
  if ( verbosity_ >= 1 ) {
    for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
      const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
      std::cout << "measuredTauLepton #" << idx << " (type = " << measuredTauLepton.type() << "): Pt = " << measuredTauLepton.pt() << ","
                << " eta = " << measuredTauLepton.eta() << " (theta = " << measuredTauLepton.p3().theta() << ")" << ", phi = " << measuredTauLepton.phi() << ","
                << " mass = " << measuredTauLepton.mass() << std::endl;
    }
  }
  double measuredMETx_rounded = roundToNdigits(measuredMETx);
  double measuredMETy_rounded = roundToNdigits(measuredMETy);
  met_.SetX(measuredMETx_rounded);
  met_.SetY(measuredMETy_rounded);
  met_.SetZ(0.0);
  TMatrixD covMET_rounded(2,2);
  covMET_rounded[0][0] = roundToNdigits(covMET[0][0]);
  covMET_rounded[1][0] = roundToNdigits(covMET[1][0]);
  covMET_rounded[0][1] = roundToNdigits(covMET[0][1]);
  covMET_rounded[1][1] = roundToNdigits(covMET[1][1]);
  if ( verbosity_ >= 1 ) {
    std::cout << "MET: Px = " << measuredMETx_rounded << ", Py = " << measuredMETy_rounded << std::endl;
    std::cout << "covMET:" << std::endl;
    covMET_rounded.Print();
    TMatrixDSym covMET_sym(2);
    covMET_sym(0,0) = covMET_rounded[0][0];
    covMET_sym(0,1) = covMET_rounded[0][1];
    covMET_sym(1,0) = covMET_rounded[1][0];
    covMET_sym(1,1) = covMET_rounded[1][1];
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

//--- determine dimension of integration space
  int idxLeg1_X = -1;
  int idxLeg1_phi = -1;
  int idxLeg1VisPtShift = -1;
  int idxLeg1_mNuNu = -1;

  int idxLeg2_X = -1;
  int idxLeg2_phi = -1;
  int idxLeg2VisPtShift = -1;
  int idxLeg2_mNuNu = -1;

  numDimensions_ = 0;

  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
    if ( idx == 0 ) {
      if ((tau1Constraint_ == true) || (tau1Constraint_ == 1)) {
         if (diTauMassConstraint_ < 0.0) {
            idxLeg1_X = numDimensions_;
            numDimensions_ += 1;
         }
      }
      else {
        idxLeg1_X = numDimensions_;
        numDimensions_ += 1;
      }
      idxLeg1_phi = numDimensions_;
      numDimensions_ += 1;
      if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay ) {
	      if ( useHadTauTF_ ) {
		    idxLeg1VisPtShift = numDimensions_;
		    ++numDimensions_;
	      }
	    } else {
	      idxLeg1_mNuNu = numDimensions_;
	      numDimensions_ += 1;
	    }
    }
    if ( idx == 1 ) {
      if ((tau1Constraint_ == false) || (tau1Constraint_ == 0)) {
       if (diTauMassConstraint_ < 0.0) {
        idxLeg2_X = numDimensions_;
        numDimensions_ += 1;
       }
      }
      else {
        idxLeg2_X = numDimensions_;
        numDimensions_ += 1;
      }
      idxLeg2_phi = numDimensions_;
      numDimensions_ += 1;
      if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay ) {
        if ( useHadTauTF_ ) {
          idxLeg2VisPtShift = numDimensions_;
          ++numDimensions_;
        }
      } else {
        idxLeg2_mNuNu = numDimensions_;
        numDimensions_ += 1;
      }
    }
  }

  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  if ( measuredTauLeptons_rounded.size() == 2 ) {
    histogramAdapter_->setMeasurement(measuredTauLeptons_rounded[0].p4(), measuredTauLeptons_rounded[1].p4(), met_);
    histogramAdapter_->bookHistograms(measuredTauLeptons_rounded[0].p4(), measuredTauLeptons_rounded[1].p4(), met_);
  }

  integrand_->setInputs(measuredTauLeptons_rounded, measuredMETx_rounded, measuredMETy_rounded, covMET_rounded);
  integrand_->setHistogramAdapter(histogramAdapter_);
#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) integrand_->enableHadTauTF();
  else integrand_->disableHadTauTF();
#endif
  integrand_->setIdxLeg1_X(idxLeg1_X);
  integrand_->setIdxLeg1_phi(idxLeg1_phi);
  integrand_->setIdxLeg1VisPtShift(idxLeg1VisPtShift);
  integrand_->setIdxLeg1_mNuNu(idxLeg1_mNuNu);
  integrand_->setIdxLeg2_X(idxLeg2_X);
  integrand_->setIdxLeg2_phi(idxLeg2_phi);
  integrand_->setIdxLeg2VisPtShift(idxLeg2VisPtShift);
  integrand_->setIdxLeg2_mNuNu(idxLeg2_mNuNu);
  integrand_->setNumDimensions(numDimensions_);
  ClassicSVfitIntegrand::gSVfitIntegrand = integrand_;

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
    treeFileName_.data());
  intAlgo_->registerCallBackFunction(*histogramAdapter_);

  //std::cout << "numDimensions = " << numDimensions_ << std::endl;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  if (idxLeg1_X != -1) {
    xl_[idxLeg1_X] = 0.;
#ifdef USE_SVFITTF
    xu_[idxLeg1_X] = 2.; // upper integration bound for x1' = visPtShift1*x1
#else
    xu_[idxLeg1_X] = 1.;
#endif
  }
  xl_[idxLeg1_phi] = 0.;
  xl_[idxLeg1_phi] = -TMath::Pi();
  xu_[idxLeg1_phi] = +TMath::Pi();
  if ( idxLeg1VisPtShift != -1 ) {
    xl_[idxLeg1VisPtShift] = 0.;
    xu_[idxLeg1VisPtShift] = 2.;
  }
  if ( idxLeg1_mNuNu != -1 ) {
    xl_[idxLeg1_mNuNu] = 0.;
    xu_[idxLeg1_mNuNu] = tauLeptonMass2;
  }
  if (idxLeg2_X != -1) {
    xl_[idxLeg2_X] = 0.;
#ifdef USE_SVFITTF
    xu_[idxLeg2_X] = 2.; // upper integration bound for x2' = visPtShift2*x2
#else
    xu_[idxLeg2_X] = 1.;
#endif
  }
  xl_[idxLeg2_phi] = -TMath::Pi();
  xu_[idxLeg2_phi] = +TMath::Pi();
  if ( idxLeg2VisPtShift != -1 ) {
    xl_[idxLeg2VisPtShift] = 0.;
    xu_[idxLeg2VisPtShift] = 2.;
  }
  if ( idxLeg2_mNuNu != -1 ) {
    xl_[idxLeg2_mNuNu] = 0.;
    xu_[idxLeg2_mNuNu] = tauLeptonMass2;
  }
  if ( verbosity_ >= 1 ) {
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << " fitParameter #" << iDimension << ": xl = " << xl_[iDimension] << ", xu = " << xu_[iDimension];
      if ( (int)iDimension == idxLeg1_X         ) std::cout << " (leg1:X)";
      if ( (int)iDimension == idxLeg1_phi       ) std::cout << " (leg1:phi)";
      if ( (int)iDimension == idxLeg1VisPtShift ) std::cout << " (leg1:VisPtShift)";
      if ( (int)iDimension == idxLeg1_mNuNu     ) std::cout << " (leg1:mNuNu)";
      if ( (int)iDimension == idxLeg2_X         ) std::cout << " (leg2:X)";
      if ( (int)iDimension == idxLeg2_phi       ) std::cout << " (leg2:phi)";
      if ( (int)iDimension == idxLeg2VisPtShift ) std::cout << " (leg2:VisPtShift)";
      if ( (int)iDimension == idxLeg2_mNuNu     ) std::cout << " (leg2:mNuNu)";
      std::cout << std::endl;
    }
  }

  double integral = 0.;
  double integralErr = 0.;
  intAlgo_->integrate(&g_C, xl_, xu_, numDimensions_, integral, integralErr);
  isValidSolution_ = histogramAdapter_->isValidSolution();

  if ( likelihoodFileName_ != "" ) {
    histogramAdapter_->writeHistograms(likelihoodFileName_);
  }

  delete [] xl_;
  delete [] xu_;

  delete intAlgo_;

  clock_->Stop("<ClassicSVfit::integrate>");
  numSeconds_cpu_ = clock_->GetCpuTime("<ClassicSVfit::integrate>");
  numSeconds_real_ = clock_->GetRealTime("<ClassicSVfit::integrate>");

  if ( verbosity_ >= 1 ) {
    clock_->Show("<ClassicSVfit::integrate>");
  }
}

void ClassicSVfit::setHistogramAdapter(classic_svFit::HistogramAdapter* histogramAdapter)
{
  if ( histogramAdapter_ ) delete histogramAdapter_;
  histogramAdapter_ = histogramAdapter;
}

classic_svFit::HistogramAdapter* ClassicSVfit::getHistogramAdapter() const
{
  return histogramAdapter_;
}
