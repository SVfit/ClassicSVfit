#ifndef TauAnalysis_ClassicSVfit_FittedTauLepton_h
#define TauAnalysis_ClassicSVfit_FittedTauLepton_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // Vector, LorentzVector

namespace classic_svFit
{
  class FittedTauLepton
  {
   public:
    /// error codes that can be read out by ClassicSVfit class
    enum ErrorCodes {
      None               = 0x00000000,
      TauDecayParameters = 0x00000001
    };

    FittedTauLepton(int iTau, int verbosity);
    ~FittedTauLepton();

    void
    setMeasuredTauLepton(const MeasuredTauLepton& measuredTauLepton);
    const MeasuredTauLepton&
    getMeasuredTauLepton() const;

    /// scale momenta of visible tau decays products
    void
    updateVisMomentum(double visPtShift);

    /// compute tau lepton momentum, given momentum of visible tau decays products and the three parameters x, nuPhi, nuMass
    void
    updateTauMomentum(double x, double phiNu, double nuMass);

    /// momentum of visible tau decay products (in labframe)  
    const LorentzVector&
    visP4() const;

    /// sum of momenta of all neutrinos produced in tau decay (in labframe)
    const LorentzVector&
    nuP4() const;

    /// tau lepton momentum (in labframe)
    const LorentzVector&
    tauP4() const;

    /// parameters x, phiNu, nuMass 
    double
    x() const;
    double
    phiNu() const;
    double
    nuMass() const;

    int
    errorCode() const;

   private:
    /// instance counter (only used for debug output)
    int iTau_;

    /// pointer to MeasuredTauLepton
    MeasuredTauLepton measuredTauLepton_;
    
    /// auxiliary data-members to speed-up numerical computations
    double measuredTauLepton_mass_;
    double measuredTauLepton_mass2_;

    /// parameters used to parametrize the tau decays
    double x_;
    double phiNu_;
    double nuMass_;

    /// error code that can be passed on
    mutable int errorCode_;

    /// momentum of visible tau decay products (in labframe)  
    mutable LorentzVector visP4_;

    /// sum of momenta of all neutrinos produced in tau decay (in labframe)
    mutable LorentzVector nuP4_;

    /// tau lepton momentum (in labframe)
    mutable LorentzVector tauP4_;

    /// local coordinate system in which momentum of visible tau decay products defines z-axis (cf. Figure 4 in arXiv:1603.05910)
    double eX_x_;
    double eX_y_;
    double eX_z_;
    double eY_x_;
    double eY_y_;
    double eY_z_;
    double eZ_x_;
    double eZ_y_;
    double eZ_z_;

    /// verbosity level      
    int verbosity_;
  };
}

#endif
