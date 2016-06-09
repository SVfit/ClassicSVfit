#ifndef TauAnalysis_ClassicSVfit_ClassicSVfit_h
#define TauAnalysis_ClassicSVfit_ClassicSVfit_h

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

class ClassicSVfit
{
 public:
  ClassicSVfit(int = 0); 
  ~ClassicSVfit();

  /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
  void addLogM_fixed(bool value, double power = 1.) 
  { 
    integrand_->addLogM_fixed(value, power); 
  }
  void addLogM_dynamic(bool value, const std::string& power = "") 
  { 
    integrand_->addLogM_dynamic(value, power); 
  }

#ifdef USE_SVFITTF
  /// set transfer functions for pT of hadronic tau decays
  void setHadTauTF(const HadTauTFBase* hadTauTF) 
  { 
    integrand_->setHadTauTF(hadTauTF);
  }
  /// enable/disable use of transfer functions for hadronic tau decays
  void enableHadTauTF() 
  { 
    integrand_->enableHadTauTF();
    useHadTauTF_ = true; 
  }
  void disableHadTauTF() 
  { 
    integrand_->disableHadTauTF();
    useHadTauTF_ = false; 
  }

  /// set correlation between hadronic tau pT and MET
  void setRhoHadTau(double rhoHadTau) 
  { 
    integrand_->setRhoHadTau(rhoHadTau);
  }
#endif

  /// number of function calls for Markov Chain integration (default is 100000)
  void setMaxObjFunctionCalls(unsigned maxObjFunctionCalls) 
  { 
    maxObjFunctionCalls_ = maxObjFunctionCalls;
  }

  /// set name of ROOT file to store histograms of di-tau pT, eta, phi, mass and transverse mass
  void setLikelihoodFileName(const std::string& likelihoodFileName)
  {
    likelihoodFileName_ = likelihoodFileName;
  }

  /// set name of ROOT file to store Markov Chain steps
  void setTreeFileName(const std::string& treeFileName)
  {
    treeFileName_ = treeFileName;
  }

  /// run integration 
  void integrate(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&);

  /// return transverse momentum (pT) of the di-tau system, uncertainty on pT of the di-tau system
  /// and maximum of likelihood function versus pT 
  double pt() const { return pt_; }
  double ptErr() const { return ptErr_; }
  double ptLmax() const { return ptLmax_; }

  /// return pseudo-rapidity (eta) of the di-tau system, uncertainty on eta of the di-tau system
  /// and maximum of likelihood function versus eta  
  double eta() const { return eta_; }
  double etaErr() const { return etaErr_; }
  double etaLmax() const { return etaLmax_; }

  /// return azimuthal angle (phi) of the di-tau system, uncertainty on phi of the di-tau system
  /// and maximum of likelihood function versus phi 
  double phi() const { return phi_; }
  double phiErr() const { return phiErr_; }
  double phiLmax() const { return phiLmax_; }

  /// return mass of the di-tau system, uncertainty on the mass of the di-tau system
  /// and maximum of likelihood function versus mass  
  double mass() const { return mass_; }
  double massErr() const { return massErr_; }
  double massLmax() const { return massLmax_; }

  /// return mass of the di-tau system, uncertainty on the mass of the di-tau system
  /// and maximum of likelihood function versus mass  
  double transverseMass() const { return transverseMass_; }
  double transverseMassErr() const { return transverseMassErr_; }
  double transverseMassLmax() const { return transverseMassLmax_; }
  
  /// return flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution() const { return isValidSolution_; }

  /// return computing time (in seconds) spent on last call to integrate method
  double getComputingTime_cpu() const { return numSeconds_cpu_; }
  double getComputingTime_real() const { return numSeconds_real_; }

 protected:

  classic_svFit::ClassicSVfitIntegrand* integrand_;
  
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;

  /// interface to Markov Chain integration algorithm
  classic_svFit::SVfitIntegratorMarkovChain* intAlgo_;
  unsigned maxObjFunctionCalls_;
  double precision_;
  std::string treeFileName_;

  /// dimension of integration region
  unsigned numDimensions_;

  /// lower and upper boundary of integration region
  double* xl_;
  double* xu_;
  
  /// histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  mutable classic_svFit::HistogramAdapter* histogramAdapter_;
  std::string likelihoodFileName_;

  /// pT, eta, phi, mass and transverse mass of di-tau system
  double pt_;
  double ptErr_;
  double ptLmax_;
  double eta_;
  double etaErr_;
  double etaLmax_;
  double phi_;
  double phiErr_;
  double phiLmax_;
  double mass_;
  double massErr_;
  double massLmax_;
  double transverseMass_;
  double transverseMassErr_;
  double transverseMassLmax_;

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
