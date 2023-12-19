#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"

#include "TauAnalysis/ClassicSVfit/interface/comp_PCA.h" // comp_PCA()

#include <TH1.h>                                         // TH1::AddDirectory()
#include <TMath.h>                                       // TMath::Nint(), TMath::Pi()

#include <cmath>                                         // atan2

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
  , useTauFlightLength_(false)
  , diTauMassConstraint_(-1.)
  , useHadTauTF_(false)
  , useStartPos_(false)
  , histogramAdapter_(nullptr)
  , intAlgo_(nullptr)
  , maxObjFunctionCalls_(100000)
  , treeFileName_("")
  , likelihoodFileName_("")
  , numDimensions_(0)
  , isValidSolution_(false)
  , clock_(nullptr)
  , numSeconds_cpu_(-1.)
  , numSeconds_real_(-1.)
  , verbosity_(verbosity)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<ClassicSVfit::ClassicSVfit>:" << std::endl;
  }

  integrand_ = new ClassicSVfitIntegrand(verbosity_);

  histogramAdapter_ = new HistogramAdapterDiTau("ditau");

  unsigned int numTaus = 2;
  legIntegrationParams_.resize(numTaus);

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
ClassicSVfit::setStartPosition(const LorentzVector& tauPlusP4, const LorentzVector& tauMinusP4)
{
  useStartPos_ = true;
  startPos_tauPlusP4_ = tauPlusP4;
  startPos_tauMinusP4_ = tauMinusP4;
}

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

  measuredEvent_ = measuredEvent;
  measuredTauLeptons_ = measuredEvent.tauLeptons();
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

  const std::vector<MeasuredMEt>& measuredMEt = measuredEvent.MEt();
  MarkovChainRecorder mcRecorder(numDimensions_);
  if ( measuredMEt.size() > 1 )
  {
    intAlgo_->registerCallBackFunction(mcRecorder);
  }

  intAlgo_->registerCallBackFunction(*histogramAdapter_);
  if ( useStartPos_ )
  {
    setStartPositionImp();
  }
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
  
  useStartPos_ = false;

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
  xl_.resize(numDimensions_);
  xh_.resize(numDimensions_);
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
  const integrationParameters& aIntParams = legIntegrationParams_[iLeg];
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

namespace
{
  double
  compStartPos_phi(const LorentzVector& startPos_tauP4, const MeasuredTauLepton& measuredTauLepton)
  {
//std::cout << "<compStartPos_phi>:" << std::endl;
    Vector beamAxis(0., 0., 1.);
    Vector eZ = normalize(measuredTauLepton.p3());
    Vector eY = normalize(compCrossProduct(eZ, beamAxis));
    Vector eX = normalize(compCrossProduct(eY, eZ));
//std::cout << "eX: x = " << eX.x() << ", y = " << eX.y() << ", z = " << eX.z() << std::endl;
//std::cout << "eY: x = " << eY.x() << ", y = " << eY.y() << ", z = " << eY.z() << std::endl;
//std::cout << "eZ: x = " << eZ.x() << ", y = " << eZ.y() << ", z = " << eZ.z() << std::endl;

    Vector nuP3 = startPos_tauP4.Vect() - measuredTauLepton.p3();
//std::cout << "nuEn = " << (startPos_tauP4 - measuredTauLepton.p4()).energy() << std::endl;
//std::cout << "nuP = " << std::sqrt(nuP3.mag2()) << std::endl;
    double nuPx_local = compScalarProduct(nuP3, eX);
    double nuPy_local = compScalarProduct(nuP3, eY);
//std::cout << "nuPx_local = " << nuPx_local << ", nuPy_local = " << nuPy_local << std::endl;

    double startPos_phi = atan2(nuPy_local, nuPx_local);
//std::cout << "startPos_phi = " << startPos_phi << std::endl; 
    return startPos_phi;
  }

  double
  compStartPos_flightLength(const MeasuredEvent& measuredEvent,
                            const LorentzVector& startPos_tauP4, const MeasuredTauLepton& measuredTauLepton)
  {
    const MeasuredHadTauDecayProduct* leadChargedHadron = measuredTauLepton.leadChargedHadron();
    assert(leadChargedHadron);

    TMatrixD decayVertexCov = measuredTauLepton.decayVertexCov() + measuredEvent.primaryVertexCov();
    bool errorFlag = false;
    TMatrixD decayVertexCovInv = invertMatrix("decayVertexCov", decayVertexCov, errorFlag);
    if ( errorFlag )
    {
      std::cerr << "ERROR: Failed to invert matrix decayVertexCov !!" << std::endl;
      assert(0);
    }

    Point pca = comp_PCA(
                  startPos_tauP4,
                  measuredTauLepton, *leadChargedHadron,
                  measuredEvent.primaryVertex(), measuredTauLepton.decayVertex(), decayVertexCovInv);
    Vector flightLength = pca - measuredEvent.primaryVertex();
    double d = std::sqrt(flightLength.mag2());

    std::pair<double,double> dmin_and_dmax = comp_dmin_and_dmax(startPos_tauP4, flightLength, decayVertexCov);
    double dmin = dmin_and_dmax.first;
    double dmax = dmin_and_dmax.second;

    double startPos_flightLength = d/(dmax - dmin);
    if ( startPos_flightLength > 1. ) startPos_flightLength = 1.;
    return startPos_flightLength;
  }

  void
  setStartPos_x(std::vector<double>& startPos_x,
                int idx, double value)
  {
    if ( idx != -1 ) startPos_x[idx] = value;
  }

  void
  setStartPos_x(std::vector<double>& startPos_x,
                const MeasuredEvent& measuredEvent,
                const LorentzVector& startPos_tauP4, const MeasuredTauLepton& measuredTauLepton,
                const integrationParameters& aIntParams)
  {
//std::cout << "<setStartPos_x>:" << std::endl;
//LorentzVector nuP4 = startPos_tauP4 - measuredTauLepton.p4();
//std::cout << "cosThetaNu = " << measuredTauLepton.p3().Dot(nuP4.Vect())/(measuredTauLepton.p()*nuP4.P()) << std::endl;
    double X = measuredTauLepton.energy()/startPos_tauP4.energy();
    double phi = compStartPos_phi(startPos_tauP4, measuredTauLepton);
    double mNuNu = (startPos_tauP4 - measuredTauLepton.p4()).mass();
    double flightLength = 0.;
    if ( aIntParams.idx_flightLength_ )
    {
      flightLength = compStartPos_flightLength(measuredEvent, startPos_tauP4, measuredTauLepton);
    }
    setStartPos_x(startPos_x, aIntParams.idx_X_,             X);
    setStartPos_x(startPos_x, aIntParams.idx_phi_,           phi);
    setStartPos_x(startPos_x, aIntParams.idx_VisPtShift_,    0.);
    setStartPos_x(startPos_x, aIntParams.idx_mNuNu_,         mNuNu);
    setStartPos_x(startPos_x, aIntParams.idx_flightLength_,  flightLength);
  }
}

void
ClassicSVfit::setStartPositionImp()
{
  int tauPlus_iLeg  = -1;
  int tauMinus_iLeg = -1;
  for ( size_t iTau = 0; iTau < measuredTauLeptons_.size(); ++iTau )
  {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[iTau];
    if      ( measuredTauLepton.charge() > 0 ) tauPlus_iLeg  = iTau;
    else if ( measuredTauLepton.charge() < 0 ) tauMinus_iLeg = iTau;
  }
  if ( tauPlus_iLeg == -1 || tauMinus_iLeg == -1 )
  {
    std::cerr << "ERROR: Failed to find decay products of tau+ and tau- !!" << std::endl;
    assert(0);
  }

  startPos_x_.resize(numDimensions_);

  const MeasuredTauLepton& measuredTauPlus  = measuredTauLeptons_[tauPlus_iLeg];
  const integrationParameters& tauPlus_aIntParams = legIntegrationParams_[tauPlus_iLeg];
  setStartPos_x(startPos_x_, measuredEvent_, startPos_tauPlusP4_, measuredTauPlus, tauPlus_aIntParams);

  const MeasuredTauLepton& measuredTauMinus  = measuredTauLeptons_[tauMinus_iLeg];
  const integrationParameters& tauMinus_aIntParams = legIntegrationParams_[tauMinus_iLeg];
  setStartPos_x(startPos_x_, measuredEvent_, startPos_tauMinusP4_, measuredTauMinus, tauMinus_aIntParams);

  // transform to match convention used by SVfitIntegratorMarkovChain::setStartPosition function
  for ( unsigned int iDimension = 0; iDimension < numDimensions_; ++iDimension )
  {
    double value = (startPos_x_[iDimension] - xl_[iDimension])/(xh_[iDimension] - xl_[iDimension]);
    if ( !(value >= 0. && value <= 1.) )
    {
      std::cerr << "WARNING: start position #" << iDimension << " outside of range [0,1] !!" << std::endl;
      if      ( value < 0. ) value = 0.;
      else if ( value > 1. ) value = 1.;
    }
    startPos_x_[iDimension] = value;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "setting start-position for Markov-Chain integration to x = { ";
    for ( unsigned int iDimension = 0; iDimension < numDimensions_; ++iDimension )
    {
      std::cout << startPos_x_[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  intAlgo_->setStartPosition(startPos_x_);
}
