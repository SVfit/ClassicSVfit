#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"

#include <TH1.h>   // TH1::AddDirectory()
#include <TMath.h> // TMath::Nint(), TMath::Pi()

using namespace classic_svFit;

namespace
{
  double g_C(const double* x, size_t dim, void* param)
  {
    return ClassicSVfitIntegrand::gSVfitIntegrand->Eval(x);
  }
}

ClassicSVfit::ClassicSVfit(int verbosity)
  : integrand_(nullptr)
  , useHadTauTF_(false)
  , useTauFlightLength_(false)
  , diTauMassConstraint_(-1.)
  , histogramAdapter_(nullptr)
  , intAlgo_(nullptr)
  , maxObjFunctionCalls_(100000)
  , treeFileName_("")
  , likelihoodFileName_("")
  , numDimensions_(0)
  , xl_(nullptr)
  , xh_(nullptr)
  , isValidSolution_(false)
  , clock_(nullptr)
  , numSeconds_cpu_(-1.)
  , numSeconds_real_(-1.)
  , verbosity_(verbosity)
{
  integrand_ = new ClassicSVfitIntegrand(verbosity_);

  histogramAdapter_ = new HistogramAdapterDiTau("ditau");

  legIntegrationParams_.resize(2);

  clock_ = new TBenchmark();
}

ClassicSVfit::~ClassicSVfit()
{
  delete integrand_;

  delete histogramAdapter_;
  for ( HistogramAdapterDiTau* histogramAdapter : histogramAdaptersMEtSystematic_ )
  {
    delete histogramAdapter;
  }
  histogramAdaptersMEtSystematic_.clear();

  if ( intAlgo_ )
  {
    delete intAlgo_;
  }

  //delete [] xl_;
  //delete [] xh_;

  delete clock_;
}

void
ClassicSVfit::enableTauFlightLength()
{
  useTauFlightLength_ = true;
  integrand_->enableTauFlightLength();
}

void
ClassicSVfit::disableTauFlightLength()
{
  useTauFlightLength_ = false;
  integrand_->disableTauFlightLength();
}

void
ClassicSVfit::enableDiTauMassConstraint(double diTauMass)
{
  diTauMassConstraint_ = diTauMass;
  integrand_->enableDiTauMassConstraint(diTauMassConstraint_);
}

void
ClassicSVfit::disableDiTauMassConstraint()
{
  diTauMassConstraint_ = -1.;
  integrand_->disableDiTauMassConstraint();
}

void
ClassicSVfit::enableLogM(double power)
{
  integrand_->enableLogM(power);
}

void
ClassicSVfit::disableLogM()
{
  integrand_->disableLogM();
}

#ifdef USE_SVFITTF
void
ClassicSVfit::enableHadTauTF(const HadTauTFBase* hadTauTF, double rhoHadTau)
{
  useHadTauTF_ = true;
  integrand_->enableHadTauTF(hadTauTF, rhoHadTau);
}

void
ClassicSVfit::disableHadTauTF()
{
  useHadTauTF_ = false;
  integrand_->disableHadTauTF();
}
#endif

void
ClassicSVfit::setHistogramAdapter(classic_svFit::HistogramAdapterDiTau* histogramAdapter)
{
  if ( histogramAdapter_ )
  {
    delete histogramAdapter_;
  }
  for ( HistogramAdapterDiTau* histogramAdapter : histogramAdaptersMEtSystematic_ )
  {
    delete histogramAdapter;
  }
  histogramAdaptersMEtSystematic_.clear();
  histogramAdapter_ = histogramAdapter;
}

const classic_svFit::HistogramAdapterDiTau*
ClassicSVfit::getHistogramAdapter(unsigned int idx) const
{
  if ( idx == 0 )
  {
    return histogramAdapter_;
  }
  else
  {
    assert(idx <= histogramAdaptersMEtSystematic_.size());
    return histogramAdaptersMEtSystematic_[idx - 1];
  }
}

void ClassicSVfit::setVerbosity(int aVerbosity)
{
  verbosity_ = aVerbosity;
  integrand_->setVerbosity(verbosity_);
}

void
ClassicSVfit::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  maxObjFunctionCalls_ = maxObjFunctionCalls;
  if ( maxObjFunctionCalls_ < 1000 )
  {
    std::cerr << "WARNING: Parameter 'maxObjFunctionCalls' = " << maxObjFunctionCalls << " too low, setting it to 1000 !!" << std::endl;
    maxObjFunctionCalls_ = 1000;
  }
}

void
ClassicSVfit::setLikelihoodFileName(const std::string& likelihoodFileName)
{
  likelihoodFileName_ = likelihoodFileName;
}

void
ClassicSVfit::setTreeFileName(const std::string& treeFileName)
{
  treeFileName_ = treeFileName;
}

void ClassicSVfit::initializeIntegrand(const MeasuredEvent& measuredEvent)
{
  integrand_->setMeasurement(measuredEvent);
  integrand_->setHistogramAdapter(histogramAdapter_);
  for ( size_t iLeg = 0; iLeg < legIntegrationParams_.size(); ++iLeg )
  {
    integrand_->initializeLegIntegrationParams(iLeg, legIntegrationParams_[iLeg]);
  }
  integrand_->setNumDimensions(numDimensions_);
  integrand_->setIntegrationRanges(xl_, xh_);
  ClassicSVfitIntegrand::gSVfitIntegrand = integrand_;
}

void ClassicSVfit::integrate(const MeasuredEvent& measuredEvent)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<ClassicSVfit::integrate>:" << std::endl;
  }

  clock_->Reset();
  clock_->Start("<ClassicSVfit::integrate>");

  measuredTauLeptons_ = measuredEvent.measuredTauLeptons();
  if ( verbosity_ >= 1 )
  {
    std::cout << measuredTauLeptons_;
    LorentzVector sumP4;
    for ( const MeasuredTauLepton& measuredTauLepton : measuredTauLeptons_ )
    {
      sumP4 += measuredTauLepton.p4();
    }
    std::cout << "visible momentum sum:" 
              << " Pt = " << sumP4.pt() << "," 
              << " phi = " << sumP4.phi() << "," 
              << " mass = " << sumP4.mass() << std::endl;
  }

  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  if ( measuredTauLeptons_.size() == 2 )
  {
    histogramAdapter_->setMeasurement(measuredEvent);
    histogramAdapter_->bookHistograms(measuredEvent);
  } else assert(0);


  initializeIntegrationParams();
  initializeIntegrand(measuredEvent);
  if ( !intAlgo_ )
  {  
    initializeIntAlgo();
  }
  intAlgo_->clearCallBackFunctions();

  const std::vector<MeasuredMEt>& measuredMEt = measuredEvent.measuredMEt();
  MarkovChainRecorder mcRecorder(numDimensions_);
  if ( measuredMEt.size() > 1 )
  {
    intAlgo_->registerCallBackFunction(mcRecorder);
  }

  intAlgo_->registerCallBackFunction(*histogramAdapter_);
  integrand_->setCentral();
  double theIntegral, theIntegralErr;
  intAlgo_->integrate(&g_C, xl_, xh_, numDimensions_, theIntegral, theIntegralErr);
  isValidSolution_ = histogramAdapter_->isValidSolution();

  if ( measuredMEt.size() > 1 )
  {
    for ( HistogramAdapterDiTau* histogramAdapter : histogramAdaptersMEtSystematic_ )
    {
      delete histogramAdapter;
    }
    histogramAdaptersMEtSystematic_.clear();
    unsigned int numMEtSystematics = measuredMEt.size() - 1;
    for ( unsigned int iMEtSystematic = 0; iMEtSystematic < numMEtSystematics; ++iMEtSystematic )
    {
      HistogramAdapterDiTau* histogramAdapter = static_cast<HistogramAdapterDiTau*>(histogramAdapter_->clone());
      integrand_->setMEtSystematic(iMEtSystematic);
      integrand_->setHistogramAdapter(histogramAdapter);

      unsigned int numPoints = mcRecorder.getNumPoints();
      for ( unsigned int iPoint = 0; iPoint < numPoints; ++iPoint )
      {
        const double* x = mcRecorder.getPoint(iPoint);
        double prob_ratio = integrand_->Eval(x)/mcRecorder.getValue(iPoint);
        histogramAdapter->fillHistograms(prob_ratio);
      }

      histogramAdaptersMEtSystematic_.push_back(histogramAdapter);
    }
  }

  if ( likelihoodFileName_ != "" )
  {
    histogramAdapter_->writeHistograms(likelihoodFileName_);
  }
  
  clock_->Stop("<ClassicSVfit::integrate>");
  numSeconds_cpu_ = clock_->GetCpuTime("<ClassicSVfit::integrate>");
  numSeconds_real_ = clock_->GetRealTime("<ClassicSVfit::integrate>");
  
  if ( verbosity_ >= 1 )
  {
    clock_->Show("<ClassicSVfit::integrate>");
  }
}

bool ClassicSVfit::isValidSolution() const 
{
  return isValidSolution_;
}

void ClassicSVfit::initializeIntegrationParams()
{
  numDimensions_ = 0;
  bool useDiTauMassConstraint = (diTauMassConstraint_ > 0);
  initializeLegIntegrationParams(0, false);
  initializeLegIntegrationParams(1, useDiTauMassConstraint);
  delete [] xl_;
  xl_ = new double[numDimensions_];
  delete [] xh_;
  xh_ = new double[numDimensions_];
  initializeLegIntegrationRanges(0);
  initializeLegIntegrationRanges(1);
  if ( verbosity_ >= 1 )
  {
    std::cout << "numDimensions = " << numDimensions_ << std::endl;
    for ( size_t iDimension = 0; iDimension < numDimensions_; ++iDimension )
    {
      std::cout << " fitParameter #" << iDimension << ": xl = " << xl_[iDimension] << ", xh = " << xh_[iDimension];
      for ( size_t iLeg = 0; iLeg < legIntegrationParams_.size(); ++iLeg )
      {
        const integrationParameters& legIntegrationParams = legIntegrationParams_[iLeg];
        if ( (int)iDimension == legIntegrationParams.idx_X_            ) std::cout << " (leg" << (iLeg + 1) << ":X)";
        if ( (int)iDimension == legIntegrationParams.idx_phi_          ) std::cout << " (leg" << (iLeg + 1) << ":phi)";
        if ( (int)iDimension == legIntegrationParams.idx_VisPtShift_   ) std::cout << " (leg" << (iLeg + 1) << ":VisPtShift)";
        if ( (int)iDimension == legIntegrationParams.idx_mNuNu_        ) std::cout << " (leg" << (iLeg + 1) << ":mNuNu)";
        if ( (int)iDimension == legIntegrationParams.idx_flightLength_ ) std::cout << " (leg" << (iLeg + 1) << ":flightLength)";
      }
      std::cout << std::endl;
    }
  }
}

double ClassicSVfit::getComputingTime_cpu() const 
{
  return numSeconds_cpu_;
}

double ClassicSVfit::getComputingTime_real() const 
{
  return numSeconds_real_;
}

void ClassicSVfit::initializeIntAlgo()
{
  //unsigned numChains = TMath::Nint(maxObjFunctionCalls_/100000.);
  unsigned numChains = 1;
  unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls_/numChains);
  unsigned numIterSampling = TMath::Nint(0.90*maxObjFunctionCalls_/numChains);
  unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.20*numIterBurnin);
  unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.60*numIterBurnin);
  if ( treeFileName_ == "" && verbosity_ >= 2 )
  {
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

void
ClassicSVfit::initializeLegIntegrationParams(size_t iLeg, bool useMassConstraint)
{
  assert(iLeg < measuredTauLeptons_.size());
  const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[iLeg];
  
  legIntegrationParams_[iLeg].reset();

  if ( measuredTauLepton.type() == MeasuredTauLepton::kPrompt )
  {
    if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay )
    {
      if ( useHadTauTF_ )
      {  
        legIntegrationParams_[iLeg].idx_VisPtShift_ = numDimensions_;
        ++numDimensions_;
      }
    }
  }
  else
  {
    if ( !useMassConstraint )
    {
      legIntegrationParams_[iLeg].idx_X_ = numDimensions_;
      ++numDimensions_;
    }
    
    legIntegrationParams_[iLeg].idx_phi_ = numDimensions_;
    ++numDimensions_;
    
    if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay )
    {
      if ( useHadTauTF_ ) 
      {
        legIntegrationParams_[iLeg].idx_VisPtShift_ = numDimensions_;
        ++numDimensions_;
      }
    }
    else
    {
      legIntegrationParams_[iLeg].idx_mNuNu_ = numDimensions_;
      ++numDimensions_;
    }

    if ( useTauFlightLength_ )
    {
      legIntegrationParams_[iLeg].idx_flightLength_ = numDimensions_;
      ++numDimensions_;
    }
  }
}

void
ClassicSVfit::initializeLegIntegrationRanges(size_t iLeg)
{
  const classic_svFit::integrationParameters& aIntParams = legIntegrationParams_[iLeg];
  if ( aIntParams.idx_X_ != -1 )
  {
    xl_[aIntParams.idx_X_] = 0.;
#ifdef USE_SVFITTF
    if ( aIntParams.idx_VisPtShift_ != -1 )
    {
      xh_[aIntParams.idx_X_] = 2.; // upper integration bound for x1' = visPtShift1*x1
    }
    else
    {
      xh_[aIntParams.idx_X_] = 1.;
    }
#else
    xh_[aIntParams.idx_X_] = 1.;
#endif
  }
  if ( aIntParams.idx_phi_ != -1 )
  {
    xl_[aIntParams.idx_phi_] = -TMath::Pi();
    xh_[aIntParams.idx_phi_] = +TMath::Pi();
  }
  if ( aIntParams.idx_VisPtShift_ != -1 )
  {
    xl_[aIntParams.idx_VisPtShift_] = 0.;
    xh_[aIntParams.idx_VisPtShift_] = 2.;
  }
  if ( aIntParams.idx_mNuNu_ != -1 )
  {
    xl_[aIntParams.idx_mNuNu_] = 0.;
    xh_[aIntParams.idx_mNuNu_] = tauLeptonMass2;
  }
  if ( aIntParams.idx_flightLength_ != -1 )
  {
    xl_[aIntParams.idx_flightLength_] = 0.;
    xh_[aIntParams.idx_flightLength_] = 1.;
  }
}
