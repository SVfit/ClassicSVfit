#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrandBase.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h" // HistogramAdapterDiTau

#include <Math/Functor.h>
#include <TMatrixD.h>

namespace classic_svFit
{
  class ClassicSVfitIntegrand : public ClassicSVfitIntegrandBase
  {
   public:
    ClassicSVfitIntegrand(int);
    ~ClassicSVfitIntegrand();

    void setDiTauMassConstraint(double diTauMass);

    /// set pointer to histograms used to keep track of pT, eta, phi, mass and transverse mass of di-tau system
    /// during Markov Chain integration
    void setHistogramAdapter(HistogramAdapterDiTau* histogramAdapter);

    /// set momenta of visible tau decay products
    void setLeptonInputs(const std::vector<classic_svFit::MeasuredTauLepton>&);

    /// evaluate Phase Space part of the integrand for given value of integration variables x
    double EvalPS(const double* x) const;

    /// evaluate the iComponent of the full integrand for given value of integration variables q.
    /// q is given in standarised range [0,1] for each dimension.
    double Eval(const double* q, unsigned int iComponent=0) const;

    /// static pointer to this (needed for interfacing the likelihood function calls to Markov Chain integration)
    static const ClassicSVfitIntegrand* gSVfitIntegrand;

   protected:
    /// momenta of visible tau decay products and of reconstructed tau leptons
    MeasuredTauLepton measuredTauLepton1_;    
    mutable FittedTauLepton fittedTauLepton1_;
    bool leg1isLeptonicTauDecay_;
    bool leg1isHadronicTauDecay_;
    bool leg1isPrompt_;
    MeasuredTauLepton measuredTauLepton2_;  
    mutable FittedTauLepton fittedTauLepton2_;
    bool leg2isLeptonicTauDecay_;
    bool leg2isHadronicTauDecay_;
    bool leg2isPrompt_;

    mutable double mVis_measured_;
    mutable double mVis2_measured_;

    double diTauMassConstraint_;
    double diTauMassConstraint2_;

    HistogramAdapterDiTau* histogramAdapter_;
  };
}

#endif
