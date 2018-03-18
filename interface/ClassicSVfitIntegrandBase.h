#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitIntegrandBase_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitIntegrandBase_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"
#endif
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

#include <Math/Functor.h>
#include <TMatrixD.h>
#include <TFormula.h>

namespace classic_svFit
{
  class ClassicSVfitIntegrandBase
  {
   public:
    /// error codes that can be read out by ClassicSVfitBase class
    enum ErrorCodes {
      None               = 0x00000000,
      MatrixInversion    = 0x00000001,
      LeptonNumber       = 0x00000010,
      TestMass           = 0x00000100,
      TauDecayParameters = 0x00001000,
    };

    ClassicSVfitIntegrandBase(int);
    virtual ~ClassicSVfitIntegrandBase();

    /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
    void addLogM_fixed(bool value, double power = 1.);
    void addLogM_dynamic(bool value, const std::string& power= "");

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
    virtual void setLeptonInputs(const std::vector<classic_svFit::MeasuredTauLepton>&);

    /// add MET  estimates, i.e. systematic effect variations
    void addMETEstimate(double, double, const TMatrixD&);

    /// remove MET estimates
    void clearMET();

    /// evaluate Phase Space part of the integrand for given value of integration variables x
    virtual double EvalPS(const double* x) const = 0;

    /// evaluate the MET TF part of the integral.
    double EvalMET_TF(double aMETx, double aMETy, const TMatrixD&) const;

    /// evaluate the MET TF part of the integral using current values of the MET variables
    /// iComponent is ans index to MET estimate, i.e. systamtic effect variation
    double EvalMET_TF(unsigned int iComponent=0) const;

    /// evaluate the iComponent of the full integrand for given value of integration variables q.
    /// q is given in standarised range [0,1] for each dimension.
    virtual double Eval(const double* q, unsigned int iComponent=0) const = 0;

    ///Transform the values fo integration variables from [0,1] to
    ///desires [xMin,xMax] range;
    void rescaleX(const double* q) const;

    int getMETComponentsSize() const;

   protected:
    /// number of tau leptons reconstructed per event
    unsigned numTaus_;

    /// momenta of reconstructed tau leptons
    std::vector<FittedTauLepton*> fittedTauLeptons_;

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

    std::vector<classic_svFit::integrationParameters> legIntegrationParams_;
    unsigned numDimensions_;
    unsigned maxNumberOfDimensions_;
    mutable double* xMin_;
    mutable double* xMax_;
    mutable double* x_;

    /// flag to enable/disable addition of log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution
    bool addLogM_fixed_;
    double addLogM_fixed_power_;
    bool addLogM_dynamic_;
    TFormula* addLogM_dynamic_formula_;

    /// error code that can be passed on
    mutable int errorCode_;

    mutable double phaseSpaceComponentCache_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
