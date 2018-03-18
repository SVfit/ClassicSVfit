#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#endif
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

#include <Math/Functor.h>
#include <TMatrixD.h>
#include <TH1.h>
#include <TString.h>
#include <TFormula.h>

namespace classic_svFit
{
  class ClassicSVfitIntegrand
  {
   public:
    /// error codes that can be read out by ClassicSVfit class
    enum ErrorCodes {
      None               = 0x00000000,
      MatrixInversion    = 0x00000001,
      LeptonNumber       = 0x00000010,
      TestMass           = 0x00000100,
      TauDecayParameters = 0x00001000,
    };

    ClassicSVfitIntegrand(int);
    ~ClassicSVfitIntegrand();

    /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
    void addLogM_fixed(bool value, double power = 1.);
    void addLogM_dynamic(bool value, const std::string& power= "");

    void setDiTauMassConstraint(double diTauMass);

    /// set pointer to histograms used to keep track of pT, eta, phi, mass and transverse mass of di-tau system
    /// during Markov Chain integration
    void setHistogramAdapter(HistogramAdapterDiTau* histogramAdapter);

    void setLegIntegrationParams(unsigned int iLeg, const classic_svFit::integrationParameters& aParams);

    void setNumDimensions(unsigned numDimensions);

    void setVerbosity(int aVerbosity);

    void setIntegrationRanges(const double* xl, const double* xu);

#ifdef USE_SVFITTF
    /// set transfer functions for pT of hadronic tau decays
    void setHadTauTF(const HadTauTFBase* hadTauTF);
    /// enable/disable use of transfer functions for hadronic tau decays
    void enableHadTauTF();
    void disableHadTauTF();

    /// set correlation between hadronic tau pT and MET
    void setRhoHadTau(double rhoHadTau);
#endif

    /// set momenta of visible tau decay products
    void setLeptonInputs(const std::vector<classic_svFit::MeasuredTauLepton>&);

    /// add MET  estimates, i.e. systematic effect variations
    void addMETEstimate(double, double, const TMatrixD&);

    /// remove MET estimates
    void clearMET();

    /// evaluate Phase Space part of the integrand for given value of integration variables x
    double EvalPS(const double* x) const;

    /// evaluate the MET TF part of the integral.
    double EvalMET_TF(double aMETx, double aMETy, const TMatrixD&) const;

    /// evaluate the MET TF part of the integral using current values of the MET variables
    /// iComponent is ans index to MET estimate, i.e. systamtic effect variation
    double EvalMET_TF(unsigned int iComponent=0) const;

    /// evaluate the iComponent of the full integrand for given value of integration variables q.
    /// q is given in standarised range [0,1] for each dimension.
    double Eval(const double* q, unsigned int iComponent=0) const;

    ///Transform the values fo integration variables from [0,1] to
    ///desires [xMin,xMax] range;
    void rescaleX(const double* q) const;

    int getMETComponentsSize() const;

    /// static pointer to this (needed for interfacing the likelihood function calls to Markov Chain integration)
    static const ClassicSVfitIntegrand* gSVfitIntegrand;

   protected:
    /// number of tau leptons reconstructed per event
    unsigned numTaus_;

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
    std::vector<FittedTauLepton*> fittedTauLeptons_;

    mutable double mVis_measured_;
    mutable double mVis2_measured_;

    /// measured MET
    std::vector<double> measuredMETx_;
    std::vector<double> measuredMETy_;

    ///MET covariance matrix
    std::vector<TMatrixD> covMET_;

    ///Inverse covariance matix elements
    double invCovMETxx_;
    double invCovMETxy_;
    double invCovMETyx_;
    double invCovMETyy_;
    double const_MET_;

#ifdef USE_SVFITTF
    /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
    std::vector<const HadTauTFBase*> hadTauTFs_;
    bool useHadTauTF_;

    double rhoHadTau_;
#endif

    classic_svFit::integrationParameters legIntegrationParams_[classic_svFit::numberOfLegs];
    unsigned numDimensions_;
    mutable double xMin_[classic_svFit::maxNumberOfDimensions];
    mutable double xMax_[classic_svFit::maxNumberOfDimensions];
    mutable double x_[classic_svFit::maxNumberOfDimensions];

    /// flag to enable/disable addition of log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution
    bool addLogM_fixed_;
    double addLogM_fixed_power_;
    bool addLogM_dynamic_;
    TFormula* addLogM_dynamic_formula_;

    double diTauMassConstraint_;
    double diTauMassConstraint2_;

    /// error code that can be passed on
    mutable int errorCode_;

    HistogramAdapterDiTau* histogramAdapter_;

    mutable double phaseSpaceComponentCache_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
