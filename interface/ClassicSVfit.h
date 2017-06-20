#ifndef TauAnalysis_ClassicSVfit_ClassicSVfit_h
#define TauAnalysis_ClassicSVfit_ClassicSVfit_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/ClassicSVfit/interface/SVfitCUBAIntegrator.h"
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#endif
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include <TBenchmark.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TMath.h>

class ClassicSVfit
{
 public:
  ClassicSVfit(int = 0);
  ~ClassicSVfit();

  /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
  void addLogM_fixed(bool value, double power = 1.);
  void addLogM_dynamic(bool value, const std::string& power = "");

  void setDiTauMassConstraint(double diTauMass);

#ifdef USE_SVFITTF
  /// set transfer functions for pT of hadronic tau decays
  void setHadTauTF(const HadTauTFBase* hadTauTF);
  /// enable/disable use of transfer functions for hadronic tau decays
  void enableHadTauTF();
  void disableHadTauTF();

  /// set correlation between hadronic tau pT and MET
  void setRhoHadTau(double rhoHadTau);
#endif

  ///use Cuba integration library instead of custom Markov Chain
  void setUseCuba(bool useCuba);

  /// number of function calls for Markov Chain integration (default is 100000)
  void setMaxObjFunctionCalls(unsigned maxObjFunctionCalls);

  /// set name of ROOT file to store histograms of di-tau pT, eta, phi, mass and transverse mass
  void setLikelihoodFileName(const std::string& likelihoodFileName);

  /// set name of ROOT file to store Markov Chain steps
  void setTreeFileName(const std::string& treeFileName);

  /// set and get histogram adapter
  void setHistogramAdapter(classic_svFit::HistogramAdapter* histogramAdapter);
  classic_svFit::HistogramAdapter* getHistogramAdapter() const;

  /// run integration
  void integrate(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&);

  /// return flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution() const;

  /// return computing time (in seconds) spent on last call to integrate method
  double getComputingTime_cpu() const;
  double getComputingTime_real() const;

  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

 protected:

   ///flag for choosing the integrator class
   bool useCuba_;

   /// initialize Matkov Chain integrator class
   void initializeMCIntegrator();

   /// print MET and its covariance matrix
   void printMET() const;

   /// print measured leptons
   void printLeptons() const;

   /// print integration range
   void printIntegrationRange() const;

    /// set integration indices and ranges for both legs
   void setIntegrationParams();

   /// set integration indices and ranges for given leg
   void setLegIntegrationParams(unsigned int iLeg);

   /// set integration ranges for given leg
   void setIntegrationRanges(unsigned int iLeg);

  classic_svFit::ClassicSVfitIntegrand* integrand_;

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
  classic_svFit::Vector met_;
  TMatrixD covMET_rounded;

  double diTauMassConstraint_ = -1.0;

  /// interface to Markov Chain integration algorithm
  classic_svFit::SVfitIntegratorMarkovChain* intAlgo_;

  /// interface to CUBA integration algorithm
  classic_svFit::SVfitCUBAIntegrator* intCubaAlgo_;

  unsigned maxObjFunctionCalls_;
  double precision_;
  std::string treeFileName_;

  /// variables indices and ranges for each leg
  classic_svFit::integrationParameters legIntegrationParams_[2];
  unsigned int numDimensions_ = 0;
  double xl_[6], xu_[6];

  /// histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  mutable classic_svFit::HistogramAdapter* histogramAdapter_;
  std::string likelihoodFileName_;

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
