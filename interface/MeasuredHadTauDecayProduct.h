#ifndef TauAnalysis_ClassicSVfitLT_MeasuredHadTauDecayProduct_h
#define TauAnalysis_ClassicSVfitLT_MeasuredHadTauDecayProduct_h

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // Vector, LorentzVector

namespace classic_svFit
{
  class MeasuredHadTauDecayProduct
  {
   public:
    MeasuredHadTauDecayProduct(int, double, double, double, double);
    MeasuredHadTauDecayProduct(const MeasuredHadTauDecayProduct&);
    ~MeasuredHadTauDecayProduct();
   
    /// return charge
    int
    charge() const;

    /// return pt in labframe
    double
    pt() const;

    /// return pseudo-rapidity in labframe
    double
    eta() const;

    /// return azimuthal angle in labframe
    double
    phi() const;

    /// return mass
    double
    mass() const;

    /// return energy in labframe
    double
    energy() const;

    /// return px in labframe
    double
    px() const;

    /// return py in labframe
    double
    py() const;
    /// return pz in labframe
    double
    pz() const;

    /// return magnitude of momentum vector in labframe
    double
    p() const;

    /// return four-vector in labframe
    LorentzVector
    p4() const;

    /// return momentum vector in labframe
    Vector
    p3() const;
  
   protected:
    /// set momentum in all coordinates systems
    void
    initialize();

   private:
    /// electric charge
    int charge_;

    /// momentum in labframe (in polar coordinates)
    double pt_;
    double eta_;
    double phi_;
    double mass_;

    /// momentum in labframe (in cartesian coordinates)
    double energy_;
    double px_;
    double py_;
    double pz_;

    /// magnitude of momentum in labframe (magnitude)
    double p_;

    /// momentum in labframe (four-vector)
    LorentzVector p4_;

    /// momentum in labframe
    Vector p3_;
  };

  // auxiliary class for sorting MeasuredHadTauDecayProduct objects
  struct sortMeasuredHadTauDecayProducts
  {
    bool
    operator() (const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct1, const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct2);
  };
}

#endif
