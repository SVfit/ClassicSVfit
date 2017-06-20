#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"

#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitCUBAIntegrator.h"

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
    double retVal = 1E16*ClassicSVfitIntegrand::gSVfitIntegrand->Eval(x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }

  int cubaIntegrand(const int *ndim, const double qq[],
                    const int *ncomp, double ff[], void *userdata)
  {

    double xx[5];
    //double xMin_[] = {0,-3.14159,0,0,-3.14159};
    //double xMax_[] = {1,3.14159,3.14159,1,3.14159};

    double xMin_[] = {0,-3.14159,0,-3.14159};
    double xMax_[] = {1,3.14159,3.14159,3.14159};


    for ( unsigned iDimension = 0; iDimension < *ndim; ++iDimension ) {
    double q_i = qq[iDimension];
    xx[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }

    ff[0] = 1E16*ClassicSVfitIntegrand::gSVfitIntegrand->Eval(xx);
    return 0;
  }
}

ClassicSVfit::ClassicSVfit(int verbosity)
  : integrand_(0),
    intAlgo_(0),
    intCubaAlgo_(0),
    maxObjFunctionCalls_(100000),
    precision_(1.e-3),
    treeFileName_(""),
    numDimensions_(0),
    histogramAdapter_(new DiTauSystemHistogramAdapter()),
    likelihoodFileName_(""),
    isValidSolution_(false),
    useHadTauTF_(false),
    useCuba_(false),
    clock_(0),
    numSeconds_cpu_(-1.),
    numSeconds_real_(-1.),
    verbosity_(verbosity)
{
  integrand_ = new ClassicSVfitIntegrand(verbosity_);
  clock_ = new TBenchmark();

  covMET_rounded.ResizeTo(2,2);

  if(!useCuba_) initializeMCIntegrator();
  else  intCubaAlgo_ = new SVfitCUBAIntegrator(verbosity_);

}

ClassicSVfit::~ClassicSVfit()
{
  delete histogramAdapter_;

  delete integrand_;

  if(intAlgo_) delete intAlgo_;
  if(intCubaAlgo_) delete intCubaAlgo_;

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

void ClassicSVfit::setUseCuba(bool useCuba) { useCuba_ = useCuba;}

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

void ClassicSVfit::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double cubaIntegral;
    double cubaIntegralErr;
    setDiTauMassConstraint(par[0]);
    intCubaAlgo_->integrate(&cubaIntegrand, xl_, xu_, numDimensions_, cubaIntegral, cubaIntegralErr);
    f = cubaIntegral;
}

void ClassicSVfit::initializeMCIntegrator(){
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
  intAlgo_->registerCallBackFunction(*histogramAdapter_);
}

void ClassicSVfit::printMET() const {

  std::cout << "MET: Px = " << met_.X() << ", Py = " << met_.Y() << std::endl;
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

void ClassicSVfit::printLeptons() const{
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
      const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
      std::cout << "measuredTauLepton #" << idx << " (type = " << measuredTauLepton.type() << "): Pt = " << measuredTauLepton.pt() << ","
                << " eta = " << measuredTauLepton.eta() << " (theta = " << measuredTauLepton.p3().theta() << ")" << ", phi = " << measuredTauLepton.phi() << ","
                << " px = "<<measuredTauLepton.px()
                << " mass = " << measuredTauLepton.mass() << std::endl;
    }
}

void ClassicSVfit::printIntegrationRange() const {

  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << " fitParameter #" << iDimension << ": xl = " << xl_[iDimension] << ", xu = " << xu_[iDimension];
      if ( (int)iDimension == legIntegrationParams_[0].idx_X_         ) std::cout << " (leg1:X)";
      if ( (int)iDimension == legIntegrationParams_[0].idx_phi_       ) std::cout << " (leg1:phi)";
      if ( (int)iDimension == legIntegrationParams_[0].idx_VisPtShift_) std::cout << " (leg1:VisPtShift)";
      if ( (int)iDimension == legIntegrationParams_[0].idx_mNuNu_     ) std::cout << " (leg1:mNuNu)";
      if ( (int)iDimension == legIntegrationParams_[1].idx_X_         ) std::cout << " (leg2:X)";
      if ( (int)iDimension == legIntegrationParams_[1].idx_phi_       ) std::cout << " (leg2:phi)";
      if ( (int)iDimension == legIntegrationParams_[1].idx_VisPtShift_) std::cout << " (leg2:VisPtShift)";
      if ( (int)iDimension == legIntegrationParams_[1].idx_mNuNu_     ) std::cout << " (leg2:mNuNu)";
      std::cout << std::endl;
    }
}

void ClassicSVfit::setIntegrationParams(){

  numDimensions_ = 0;
  legIntegrationParams_[0].reset();
   legIntegrationParams_[1].reset();
   setLegIntegrationParams(0);
   setLegIntegrationParams(1);
}

void ClassicSVfit::setLegIntegrationParams(unsigned int iLeg){

    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[iLeg];

    if (!(iLeg == 1  && diTauMassConstraint_ > 0.0)) legIntegrationParams_[iLeg].idx_X_ = numDimensions_++;
    legIntegrationParams_[iLeg].idx_phi_ = numDimensions_++;
    if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay) {
	      if ( useHadTauTF_ ) legIntegrationParams_[iLeg].idx_VisPtShift_ = numDimensions_++;
	  } else {
	      legIntegrationParams_[iLeg].idx_mNuNu_ = numDimensions_++;
	    }
      setIntegrationRanges(iLeg);
}

void ClassicSVfit::setIntegrationRanges(unsigned int iLeg){

  const classic_svFit::integrationParameters & aIntParams = legIntegrationParams_[iLeg];

  xl_[aIntParams.idx_X_] = 0.;
#ifdef USE_SVFITTF
  xu_[aIntParams.idx_X_] = 2.; // upper integration bound for x1' = visPtShift1*x1
#else
  xu_[aIntParams.idx_X_] = 1.;
#endif
  xl_[aIntParams.idx_phi_] = -TMath::Pi();
  xu_[aIntParams.idx_phi_] = +TMath::Pi();
  if ( aIntParams.idx_VisPtShift_ != -1 ) {
    xl_[aIntParams.idx_VisPtShift_] = 0.;
    xu_[aIntParams.idx_VisPtShift_] = 2.;
  }
  if ( aIntParams.idx_mNuNu_ != -1 ) {
    xl_[aIntParams.idx_mNuNu_] = 0.;
    xu_[aIntParams.idx_mNuNu_] = tauLeptonMass2;
  }
}

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
ClassicSVfit::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons,
                        double measuredMETx, double measuredMETy,
                        const TMatrixD& covMET)
{
  if ( verbosity_ >= 1 ) std::cout << "<ClassicSVfit::integrate>:" << std::endl;

  clock_->Reset();
  clock_->Start("<ClassicSVfit::integrate>");

  measuredTauLeptons_ = measuredTauLeptons;
  for (std::vector<MeasuredTauLepton>::iterator measuredTauLepton = measuredTauLeptons_.begin();
  measuredTauLepton != measuredTauLeptons_.end(); ++measuredTauLepton ) measuredTauLepton->roundToNdigits();
  std::sort(measuredTauLeptons_.begin(), measuredTauLeptons_.end(), sortMeasuredTauLeptons());
  if ( verbosity_ >= 1 ) printLeptons();

  met_.SetX(roundToNdigits(measuredMETx));
  met_.SetY(roundToNdigits(measuredMETy));
  met_.SetZ(0.0);
  covMET_rounded[0][0] = roundToNdigits(covMET[0][0]);
  covMET_rounded[1][0] = roundToNdigits(covMET[1][0]);
  covMET_rounded[0][1] = roundToNdigits(covMET[0][1]);
  covMET_rounded[1][1] = roundToNdigits(covMET[1][1]);
  if ( verbosity_ >= 1 ) printMET();

  setIntegrationParams();
  if ( verbosity_ >= 1 ) printIntegrationRange();

  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  if ( measuredTauLeptons_.size() == 2 ) {
    histogramAdapter_->setMeasurement(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), met_);
    histogramAdapter_->bookHistograms(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), met_);
  }

  integrand_->setInputs(measuredTauLeptons_, met_.X(), met_.Y(), covMET_rounded);
  integrand_->setHistogramAdapter(histogramAdapter_);
#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) integrand_->enableHadTauTF();
  else integrand_->disableHadTauTF();
#endif
  integrand_->setLegIntegrationParams(0,legIntegrationParams_[0]);
  integrand_->setLegIntegrationParams(1,legIntegrationParams_[1]);
  integrand_->setNumDimensions(numDimensions_);
  ClassicSVfitIntegrand::gSVfitIntegrand = integrand_;

  double integral = 0.;
  double integralErr = 0.;
  intAlgo_->integrate(&g_C, xl_, xu_, numDimensions_, integral, integralErr);

/*
/////TEST
TFile file("CubaTest.root","RECREATE");
double visMass = ( measuredTauLeptons_[0].p4() + measuredTauLeptons_[1].p4()).mass();
double minMass = visMass/1.0125;
double maxMass = TMath::Max(1.e+4, 1.e+1*minMass);
TH1D *h = (TH1D*)HistogramTools::makeHistogram("Cuba_histogramMass", minMass, maxMass, 1.025);

  double cubaIntegral;
  double cubaIntegralErr;

  maxMass = 0;
  float maxIntegral = 0;
  for(unsigned int iMassPoint=1;iMassPoint<h->GetNbinsX();++iMassPoint){
    float testMass = h->GetBinCenter(iMassPoint);
    setDiTauMassConstraint(testMass);
    intCubaAlgo_->integrate(&cubaIntegrand, xl_, xu_, numDimensions_, cubaIntegral, cubaIntegralErr);

    h->SetBinContent(iMassPoint, cubaIntegral);
    if(cubaIntegral>maxIntegral){
      maxIntegral = cubaIntegral;
      maxMass = testMass;
    }
  }
  file.Write();
  std::cout<<"CUBA integral. Mass for max integral: "
            <<maxMass<<" integral: "<<maxIntegral<<std::endl;
            */

  if ( likelihoodFileName_ != "" ) {
    histogramAdapter_->writeHistograms(likelihoodFileName_);
  }

  clock_->Stop("<ClassicSVfit::integrate>");
  numSeconds_cpu_ = clock_->GetCpuTime("<ClassicSVfit::integrate>");
  numSeconds_real_ = clock_->GetRealTime("<ClassicSVfit::integrate>");

  if ( verbosity_ >= 1 ) {
    clock_->Show("<ClassicSVfit::integrate>");
  }
}

void ClassicSVfit::setHistogramAdapter(classic_svFit::HistogramAdapter* histogramAdapter)
{
  if (histogramAdapter_) delete histogramAdapter_;
  histogramAdapter_ = histogramAdapter;
}

classic_svFit::HistogramAdapter* ClassicSVfit::getHistogramAdapter() const
{
  return histogramAdapter_;
}
