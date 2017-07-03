#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
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
      None            = 0x00000000,
      MatrixInversion = 0x00000001,
      LeptonNumber    = 0x00000010,
      TestMass        = 0x00000100
    };

    ClassicSVfitIntegrand(int);
    ~ClassicSVfitIntegrand();

    /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
    void addLogM_fixed(bool value, double power = 1.);
    void addLogM_dynamic(bool value, const std::string& power= "");

    void setDiTauMassConstraint(double diTauMass);

    /// set pointer to histograms used to keep track of pT, eta, phi, mass and transverse mass of di-tau system
    /// during Markov Chain integration
    void setHistogramAdapter(HistogramAdapter* histogramAdapter);

    void setLegIntegrationParams(unsigned int iLeg,
                                 const classic_svFit::integrationParameters & aParams);

    void setNumDimensions(unsigned numDimensions);

    void setIntegrationRanges(const double* xl, const double* xu);

    void computeVisMom(const double & visPtShift1, const double & visPtShift2);

#ifdef USE_SVFITTF
    /// set transfer functions for pT of hadronic tau decays
    void setHadTauTF(const HadTauTFBase* hadTauTF);
    /// enable/disable use of transfer functions for hadronic tau decays
    void enableHadTauTF();
    void disableHadTauTF();

    /// set correlation between hadronic tau pT and MET
    void setRhoHadTau(double rhoHadTau);
#endif

    /// set momenta of visible tau decay products and of reconstructed missing transverse energy
    void setInputs(const std::vector<classic_svFit::MeasuredTauLepton>&, double, double, const TMatrixD&);

    /// evaluate Phase Space part of the integrand for given value of integration variables x
    double EvalPS(const double* x) const;

    /// evaluate the MET TF part of the integral.
    double EvalMET_TF(const double & aMETx, const double & aMETy, const TMatrixD&) const;

    /// evaluate the MET TF part of the integral using current values of the MET variables
    double EvalMET_TF() const;

    /// evaluate the full integrand for given value of integration variables q.
    /// q is given in standarised range [0,1] for each dimension.
    double Eval(const double* q) const;

    ///Transform the values fo integration variables from [0,1] to
    ///desires [xMin,xMax] range;
    void rescaleX(const double* q) const;

    /// static pointer to this (needed for interfacing the likelihood function calls to Markov Chain integration)
    static const ClassicSVfitIntegrand* gSVfitIntegrand;

   protected:

    mutable LorentzVector vis1P4_, vis2P4_;
    mutable LorentzVector nu1P4_, nu2P4_;
    mutable LorentzVector tau1P4_, tau2P4_;
    mutable double vis1P_, vis2P_;
    mutable double vis1En_, vis2En_;

    /// measured tau leptons
    MeasuredTauLepton measuredTauLepton1_;
    bool leg1isLep_;
    double leg1Mass_;
    double leg1Mass2_;
    double leg1eX_x_;
    double leg1eX_y_;
    double leg1eX_z_;
    double leg1eY_x_;
    double leg1eY_y_;
    double leg1eY_z_;
    double leg1eZ_x_;
    double leg1eZ_y_;
    double leg1eZ_z_;
    MeasuredTauLepton measuredTauLepton2_;
    bool leg2isLep_;
    double leg2Mass_;
    double leg2Mass2_;
    double leg2eX_x_;
    double leg2eX_y_;
    double leg2eX_z_;
    double leg2eY_x_;
    double leg2eY_y_;
    double leg2eY_z_;
    double leg2eZ_x_;
    double leg2eZ_y_;
    double leg2eZ_z_;

    mutable double mVis_measured_;
    mutable double mVis2_measured_;

    Vector beamAxis_;

    /// measured MET
    double measuredMETx_;
    double measuredMETy_;

    /// inverse of MET covariance matrix
    TMatrixD invCovMET_;
    double invCovMETxx_;
    double invCovMETxy_;
    double invCovMETyx_;
    double invCovMETyy_;
    double const_MET_;

#ifdef USE_SVFITTF
    /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
    const HadTauTFBase* hadTauTF1_;
    const HadTauTFBase* hadTauTF2_;
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

    double diTauMassConstraint_ = -1.0;

    /// error code that can be passed on
    int errorCode_;

    HistogramAdapter* histogramAdapter_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
