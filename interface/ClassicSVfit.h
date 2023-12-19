#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitBase_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitBase_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"      // ClassicSVfitIntegrand
#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"              // MeasuredEvent
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"          // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"      // HistogramAdapterDiTau*
#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h" // SVfitIntegratorMarkovChain
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"                    // HadTauTFBase
#endif

#include <TBenchmark.h>                                                    // TBenchmark

class ClassicSVfit
{
 public:
  ClassicSVfit(int = 0);
  ~ClassicSVfit();

  /// enable/disable use of transfer functions for (default is disabled)
  void
  enableTauFlightLength();
  void
  disableTauFlightLength();

  /// contrain tau-pair mass to given value
  void
  enableDiTauMassConstraint(double diTauMass);
  void
  disableDiTauMassConstraint();

  /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is disabled)
  void
  enableLogM(double power = 1.);
  void
  disableLogM();

#ifdef USE_SVFITTF
  /// enable/disable use of transfer functions for pT of hadronic tau decays;
  /// the parameter rhoHadTau specifies the correlation between hadronic tau pT and MET
  void
  enableHadTauTF(const HadTauTFBase* hadTauTF, double rhoHadTau = 0.);
  void
  disableHadTauTF();
#endif

  /// set start-position for Markov-Chain integration;
  /// calling this function is not neccessary, but reduces the computing time 
  /// by avoiding that SVfitIntegratorMarkovChain needs to find the start-position by try & error
  void
  setStartPosition(const classic_svFit::LorentzVector& tauPlusP4, const classic_svFit::LorentzVector& tauMinusP4);

  /// set and get histogram adapter
  void
  setHistogramAdapter(classic_svFit::HistogramAdapterDiTau* histogramAdapter);
  const classic_svFit::HistogramAdapterDiTau*
  getHistogramAdapter(unsigned int idx = 0) const;

  ///set verbosity level.
  ///Level 0 - mute, level 1 - print inputs, level 2 - print integration details
  void
  setVerbosity(int aVerbosity);

  /// number of function calls for Markov Chain integration (default is 100000)
  void
  setMaxObjFunctionCalls(unsigned maxObjFunctionCalls);

  /// set name of ROOT file to store histograms of di-tau pT, eta, phi, mass and transverse mass
  void
  setLikelihoodFileName(const std::string& likelihoodFileName);

  /// set name of ROOT file to store Markov Chain steps
  void
  setTreeFileName(const std::string& treeFileName);

  /// prepare the integrand
  void
  initializeIntegrand(const classic_svFit::MeasuredEvent& measuredEvent);

  /// run integration with Markov Chain
  void
  integrate(const classic_svFit::MeasuredEvent& measuredEvent);

  /// return flag indicating if algorithm succeeded to find valid solution
  bool
  isValidSolution() const;

  /// return computing time (in seconds) spent on last call to integrate method
  double
  getComputingTime_cpu() const;
  double
  getComputingTime_real() const;

 protected:
  /// initialize Markov Chain integrator class
  void
  initializeIntAlgo();

  /// initialize integration indices and ranges 
  void
  initializeIntegrationParams();
  void
  initializeLegIntegrationParams(size_t iLeg, bool useDiTauMassConstraint);
  void
  initializeLegIntegrationRanges(size_t iLeg);

  /// compute start-position for Markov-Chain integration
  /// based on the tau+ and tau- four-vectors passed to setStartPosition function
  void
  setStartPositionImp();

  classic_svFit::ClassicSVfitIntegrand* integrand_;

  /// reference to MeasuredTauLepton objects of MeasuredEvent given as parameter to integrate() function
  classic_svFit::MeasuredEvent measuredEvent_;
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;

  /// enable use of transverse impact parameter (for 1-prongs) and tau decay vertex (for 3-prongs)
  bool useTauFlightLength_;

  /// tau-pair mass contraint (default is -1, corresponding to constraint being disabled)
  double diTauMassConstraint_;

  /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
  bool useHadTauTF_; 

  /// set start-position for Markov-Chain integration
  bool useStartPos_;
  classic_svFit::LorentzVector startPos_tauPlusP4_;
  classic_svFit::LorentzVector startPos_tauMinusP4_;
  std::vector<double> startPos_x_;

  /// histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  mutable classic_svFit::HistogramAdapterDiTau* histogramAdapter_;
  mutable std::vector<classic_svFit::HistogramAdapterDiTau*> histogramAdaptersMEtSystematic_;

  /// interface to Markov Chain integration algorithm
  classic_svFit::SVfitIntegratorMarkovChain* intAlgo_;
  unsigned int maxObjFunctionCalls_;
  std::string treeFileName_;
  std::string likelihoodFileName_;

  /// variables indices and ranges for each leg
  std::vector<classic_svFit::integrationParameters> legIntegrationParams_;
  unsigned int numDimensions_;
  std::vector<double> xl_;
  std::vector<double> xh_;

  /// flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;
  double numSeconds_cpu_;
  double numSeconds_real_;

  /// verbosity level
  int verbosity_;
};

#endif
