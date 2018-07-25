#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
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
  double g_C(const double* x, size_t dim, void* param)
  {
    return ClassicSVfitIntegrand::gSVfitIntegrand->Eval(x);
  }
}

ClassicSVfit::ClassicSVfit(int verbosity)
  : ClassicSVfitBase(verbosity)
  , diTauMassConstraint_(-1.)
  , histogramAdapter_(new HistogramAdapterDiTau("ditau"))
{
  integrand_ = new ClassicSVfitIntegrand(verbosity_);
  legIntegrationParams_.resize(2);
  xl_ = new double[6];
  xh_ = new double[6];
}

ClassicSVfit::~ClassicSVfit()
{}

void ClassicSVfit::setDiTauMassConstraint(double diTauMass)
{
  diTauMassConstraint_ = diTauMass;
  (static_cast<ClassicSVfitIntegrand*>(integrand_))->setDiTauMassConstraint(diTauMassConstraint_);
}

void ClassicSVfit::initializeMCIntegrator()
{
  ClassicSVfitBase::initializeMCIntegrator();
  intAlgo_->registerCallBackFunction(*histogramAdapter_);
}

void ClassicSVfit::setIntegrationParams(bool useDiTauMassConstraint)
{
  numDimensions_ = 0;
  legIntegrationParams_[0].reset();
  legIntegrationParams_[1].reset();
  setLegIntegrationParams(0, false);
  setLegIntegrationParams(1, useDiTauMassConstraint);
  if ( verbosity_ >= 1 ) printIntegrationRange();
}

void ClassicSVfit::prepareIntegrand()
{
  integrand_->setLeptonInputs(measuredTauLeptons_);
  (static_cast<ClassicSVfitIntegrand*>(integrand_))->setHistogramAdapter(histogramAdapter_);
#ifdef USE_SVFITTF
  if ( useHadTauTF_ ) integrand_->enableHadTauTF();
  else integrand_->disableHadTauTF();
#endif
  for ( unsigned iLeg = 0; iLeg < legIntegrationParams_.size(); ++iLeg ) {
    integrand_->setLegIntegrationParams(iLeg, legIntegrationParams_[iLeg]);
  }
  integrand_->setNumDimensions(numDimensions_);
  integrand_->setIntegrationRanges(xl_, xh_);
  ClassicSVfitIntegrand::gSVfitIntegrand = static_cast<ClassicSVfitIntegrand*>(integrand_);
}

void ClassicSVfit::prepareLeptonInput(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  measuredTauLeptons_ = measuredTauLeptons;
  for (std::vector<MeasuredTauLepton>::iterator measuredTauLepton = measuredTauLeptons_.begin();
       measuredTauLepton != measuredTauLeptons_.end(); ++measuredTauLepton ) measuredTauLepton->roundToNdigits();
  std::sort(measuredTauLeptons_.begin(), measuredTauLeptons_.end(), sortMeasuredTauLeptons());
  if ( verbosity_ >= 1 ) {
    printLeptons();
    LorentzVector sumP4;
    for (std::vector<MeasuredTauLepton>::iterator measuredTauLepton = measuredTauLeptons_.begin();
	 measuredTauLepton != measuredTauLeptons_.end(); ++measuredTauLepton ) {
      sumP4 += measuredTauLepton->p4();
    }
    std::cout << "visible momentum sum: Pt = " << sumP4.pt() << ", phi = " << sumP4.phi() << ", mass = " << sumP4.mass() << std::endl;
  }
}

void ClassicSVfit::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons,
			     double measuredMETx, double measuredMETy,
			     const TMatrixD& covMET)
{
  if ( verbosity_ >= 1 ) std::cout << "<ClassicSVfit::integrate>:" << std::endl;

  clock_->Reset();
  clock_->Start("<ClassicSVfit::integrate>");

  prepareLeptonInput(measuredTauLeptons);
  integrand_->clearMET();
  addMETEstimate(measuredMETx, measuredMETy, covMET);
  bool useDiTauMassConstraint = (diTauMassConstraint_ > 0);
  setIntegrationParams(useDiTauMassConstraint);
  prepareIntegrand();
  if ( !intAlgo_ ) initializeMCIntegrator();

  // CV: book histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  if ( measuredTauLeptons_.size() == 2 ) {
    met_.SetX(measuredMETx);
    met_.SetY(measuredMETy);
    histogramAdapter_->setMeasurement(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), met_);
    histogramAdapter_->bookHistograms(measuredTauLeptons_[0].p4(), measuredTauLeptons_[1].p4(), met_);
  } else assert(0);
  
  double theIntegral, theIntegralErr;
  intAlgo_->integrate(&g_C, xl_, xh_, numDimensions_, theIntegral, theIntegralErr);
  isValidSolution_ = histogramAdapter_->isValidSolution();
  
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

void ClassicSVfit::setHistogramAdapter(classic_svFit::HistogramAdapterDiTau* histogramAdapter)
{
  if ( histogramAdapter_ ) delete histogramAdapter_;
  histogramAdapter_ = histogramAdapter;
}

classic_svFit::HistogramAdapterDiTau* ClassicSVfit::getHistogramAdapter() const
{
  return histogramAdapter_;
}
