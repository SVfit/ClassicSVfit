#ifndef TauAnalysis_ClassicSVfit_ClassicSVfit_h
#define TauAnalysis_ClassicSVfit_ClassicSVfit_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitBase.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

class ClassicSVfit : public ClassicSVfitBase
{
 public:
  ClassicSVfit(int = 0);
  ~ClassicSVfit();

  void setDiTauMassConstraint(double diTauMass);

  /// set and get histogram adapter
  void setHistogramAdapter(classic_svFit::HistogramAdapterDiTau* histogramAdapter);
  classic_svFit::HistogramAdapterDiTau* getHistogramAdapter() const;

  /// prepare the integrand
  void prepareIntegrand();

  /// run integration with Markov Chain
  void integrate(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&);

 protected:
  /// initialize Markov Chain integrator class
  void initializeMCIntegrator();

  /// set integration indices and ranges for both legs
  /// when useMassConstraint is true reduce number of
  /// dimension by using the mass contraint
  void setIntegrationParams(bool useMassConstraint=false);

  double diTauMassConstraint_;

  /// histograms for evaluation of pT, eta, phi, mass and transverse mass of di-tau system
  mutable classic_svFit::HistogramAdapterDiTau* histogramAdapter_;
};

#endif
