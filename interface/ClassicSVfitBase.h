#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitBase_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitBase_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#endif
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TBenchmark.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TMath.h>

class ClassicSVfitBase
{
 public:
  ClassicSVfitBase(int = 0);
  virtual ~ClassicSVfitBase();

  /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
  void addLogM_fixed(bool value, double power = 1.);
  void addLogM_dynamic(bool value, const std::string& power = "");

#ifdef USE_SVFITTF
  /// set transfer functions for pT of hadronic tau decays
  void setHadTauTF(const HadTauTFBase* hadTauTF);
  /// enable/disable use of transfer functions for hadronic tau decays
  void enableHadTauTF();
  void disableHadTauTF();

  /// set correlation between hadronic tau pT and MET
  void setRhoHadTau(double rhoHadTau);
#endif

  ///set verbosity level.
  ///Level 0 - mute, level 1 - print inputs, level 2 - print integration details
  void setVerbosity(int aVerbosity);

  /// number of function calls for Markov Chain integration (default is 100000)
  void setMaxObjFunctionCalls(unsigned maxObjFunctionCalls);

  /// set name of ROOT file to store histograms of di-tau pT, eta, phi, mass and transverse mass
  void setLikelihoodFileName(const std::string& likelihoodFileName);

  /// set name of ROOT file to store Markov Chain steps
  void setTreeFileName(const std::string& treeFileName);

  /// prepare the integrand
  virtual void prepareIntegrand() = 0;

  /// prepare input measurements
  virtual void prepareLeptonInput(const std::vector<classic_svFit::MeasuredTauLepton>& measuredTauLeptons) = 0;

  /// add MET estimate, i.e systematic effect variation
  void addMETEstimate(double measuredMETx, double measuredMETy, const TMatrixD& covMET);

  /// remove MET estimated from integrand
  void clearMET();

  /// run integration with Markov Chain
  virtual void integrate(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&) = 0;

  /// return maximum of integrand within integration domain
  double getProbMax() const { return intAlgo_->getProbMax(); }

  /// return flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution() const;

  /// return computing time (in seconds) spent on last call to integrate method
  double getComputingTime_cpu() const;
  double getComputingTime_real() const;

 protected:
  ///flag for choosing the integrator class
  bool useCuba_;

  /// initialize Markov Chain integrator class
  virtual void initializeMCIntegrator();

  /// print MET and its covariance matrix
  void printMET(double measuredMETx, double measuredMETy, const TMatrixD& covMET) const;

  /// print measured leptons
  void printLeptons() const;

  /// print integration range
  void printIntegrationRange() const;

  /// set integration indices and ranges for given leg
  void setLegIntegrationParams(unsigned iLeg, bool useMassConstraint=false);
  
  /// set integration ranges for given leg
  void setIntegrationRanges(unsigned iLeg);

  classic_svFit::ClassicSVfitIntegrandBase* integrand_;

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
  classic_svFit::Vector met_;

  /// interface to Markov Chain integration algorithm
  classic_svFit::SVfitIntegratorMarkovChain* intAlgo_;
  unsigned maxObjFunctionCalls_;
  std::string treeFileName_;
  std::string likelihoodFileName_;

  /// variables indices and ranges for each leg
  std::vector<classic_svFit::integrationParameters> legIntegrationParams_;
  unsigned numDimensions_;
  double* xl_;
  double* xh_;

  /// flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution_;

  /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
  bool useHadTauTF_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;
  double numSeconds_cpu_;
  double numSeconds_real_;

  /// verbosity level
  int verbosity_;
};

#endif
