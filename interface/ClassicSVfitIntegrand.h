#ifndef TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h
#define TauAnalysis_ClassicSVfit_ClassicSVfitIntegrand_h

#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"       // FittedTauLepton
#ifdef USE_SVFITTF
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"               // HadTauTFBase
#endif
#include "TauAnalysis/ClassicSVfit/interface/MarkovChainRecorder.h"   // MarkovChainRecorder
#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"         // MeasuredEvent 
#include "TauAnalysis/ClassicSVfit/interface/MeasuredMEt.h"           // MeasuredMEt
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"     // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"     // integrationParameters
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h" // HistogramAdapterDiTau

#include <Math/Functor.h>
#include <TMatrixD.h>                                                 // TMatrixD
#include <TFormula.h>                   

namespace classic_svFit
{
  class ClassicSVfitIntegrand
  {
   public:
    /// error codes that can be read out by ClassicSVfit class
    enum ErrorCodes
    {
      None                     = 0x00000000,
      MatrixInversion          = 0x00000001,
      LeptonNumber             = 0x00000010,
      MissingVertex            = 0x00000100,
      MissingLeadChargedHadron = 0x00001000,
      TauDecayParameters       = 0x00010000
    };

    ClassicSVfitIntegrand(int);
    ~ClassicSVfitIntegrand();

    /// switch between computation of "central" values and computation of systematic uncertainties;
    /// in case of "central" values, all elements of the likelihood function are evaluated,
    /// while for systematic uncertainties only the MEtTF is evaluated to reduce computing time.
    void
    setCentral();
    void
    setMEtSystematic(unsigned int idx);

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

    /// set pointer to histograms used to keep track of pT, eta, phi, mass and transverse mass of di-tau system
    /// during Markov Chain integration
    void
    setHistogramAdapter(HistogramAdapterDiTau* histogramAdapter);

    void
    initializeLegIntegrationParams(size_t iLeg, const integrationParameters& aParams);

    void
    setNumDimensions(unsigned int numDimensions);

    void
    setVerbosity(int aVerbosity);

    void
    setIntegrationRanges(const double* xl, const double* xh);

#ifdef USE_SVFITTF
    /// enable/disable use of transfer functions for pT of hadronic tau decays;
    /// the parameter rhoHadTau specifies the correlation between hadronic tau pT and MET
    void
    enableHadTauTF(const HadTauTFBase* hadTauTF, double rhoHadTau = 0.);
    void
    disableHadTauTF();
#endif

    /// set momenta of visible tau decay products
    void
    setMeasurement(const MeasuredEvent& measuredEvent);

    /// evaluate integrand for given value of integration variables q,
    /// where q is given in standardized range [0,1] for each dimension
    double
    Eval(const double* q) const;

    /// static pointer to this (needed for interfacing the likelihood function calls to Markov Chain integration)
    static const ClassicSVfitIntegrand* gSVfitIntegrand;

   protected:
    /// transform the values fo integration variables from [0,1] to desired [xMin,xMax] range
    void
    rescaleX(const double* q) const;

    /// evaluate phase-space part of the integrand for given value of integration variables x
    double
    EvalPS() const;

    /// evaluate part of the integrand related to transverse impact parameter (for 1-prongs) and tau decay vertex (for 3-prongs) 
    /// and to exponential decay probability of tau lepton as function of flight length
    double
    EvalFlightLength() const;

    /// evaluate the MET TF part of the integral.
    double
    EvalMEtTF(const MeasuredMEt& measuredMEt) const;

    /// number of tau leptons reconstructed per event
    unsigned int numTaus_;

    /// momenta of reconstructed tau leptons
    mutable FittedTauLepton fittedTauLepton1_;
    mutable FittedTauLepton fittedTauLepton2_;
    std::vector<FittedTauLepton*> fittedTauLeptons_;

    /// momenta and mass of visible tau decay products
    MeasuredEvent measuredEvent_;
    Point primaryVertex_;
    MeasuredTauLepton measuredTauLepton1_;
    Point leg1decayVertex_;
    TMatrixD leg1decayVertexCov_;
    TMatrixD leg1decayVertexCovInv_;
    const MeasuredHadTauDecayProduct* leg1leadChargedHadron_;
    bool leg1isLeptonicTauDecay_;
    bool leg1isHadronicTauDecay_;
    bool leg1isPrompt_;
    MeasuredTauLepton measuredTauLepton2_;
    Point leg2decayVertex_;
    TMatrixD leg2decayVertexCov_;
    TMatrixD leg2decayVertexCovInv_;
    const MeasuredHadTauDecayProduct* leg2leadChargedHadron_;
    bool leg2isLeptonicTauDecay_;
    bool leg2isHadronicTauDecay_;
    bool leg2isPrompt_;

    double mVis_measured_;
    double mVis2_measured_;

    /// flags to switch between computation of "central" values and computation of systematic uncertainties;
    /// in case of "central" values, all elements of the likelihood function are evaluated,
    /// while for systematic uncertainties only the MEtTF is evaluated to reduce computing time.
    bool isCentral_;
    unsigned int idxMEtSystematic_;

    /// flag to enable/disable use of transverse impact parameter (for 1-prongs) and tau decay vertex (for 3-prongs) information
    bool useTauFlightLength_;
    double const_FlightLength1_;
    double const_FlightLength2_;

    /// flag to enable/disable tau-pair mass constraint
    double diTauMassConstraint_;
    double diTauMassConstraint2_;

    /// flag to enable/disable use of log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution
    bool addLogM_;
    double addLogM_power_;

#ifdef USE_SVFITTF
    /// account for resolution on pT of hadronic tau decays via appropriate transfer functions;
    /// the parameter rhoHadTau specifies the correlation between hadronic tau pT and MET
    std::vector<const HadTauTFBase*> hadTauTFs_;
    double rhoHadTau_;
    bool useHadTauTF_;
#endif

    std::vector<classic_svFit::integrationParameters> legIntegrationParams_;
    unsigned int numDimensions_;
    unsigned int maxNumberOfDimensions_;
    mutable double* xMin_;
    mutable double* xMax_;
    mutable double* x_;

    /// error code that can be passed on
    mutable unsigned int errorCode_;

    /// cached return values of EvalPS() and EvalFlightLength() functions;
    /// these functions are called for the "central" value only (not for systematic uncertainties)
    /// to reduce the computing time
    mutable double probPS_;
    mutable double probFlightLength_;

    HistogramAdapterDiTau* histogramAdapter_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
