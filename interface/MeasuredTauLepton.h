#ifndef TauAnalysis_ClassicSVfit_MeasuredTauLepton_h
#define TauAnalysis_ClassicSVfit_MeasuredTauLepton_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredHadTauDecayProduct.h" // MeasuredHadTauDecayProduct
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"          // Vector, LorentzVector, Point

#include <string>                                                          // std::string
#include <ostream>                                                         // std::ostream
#include <vector>                                                          // std::vector<>

namespace classic_svFit
{
  class MeasuredTauLepton
  {
   public:
    /**
       \enum    MeasuredTauLepton::kDecayType
       \brief   enumeration of all tau decay types
    */
    enum kDecayType
    {
      kUndefinedDecayType,
      kTauToHadDecay,  /* < hadronic tau lepton decay                                        */
      kTauToElecDecay, /* < tau lepton decay to electron                                     */
      kTauToMuDecay,   /* < tau lepton decay to muon                                         */
      kPrompt          /* < electron or muon directly originating from LFV Higgs boson decay */
    };

    MeasuredTauLepton();
    MeasuredTauLepton(int,
                      int, double, double, double, double,
                      int = -1, const std::vector<MeasuredHadTauDecayProduct>* = 0);
    MeasuredTauLepton(int, 
                      int, double, double, double, double, 
                      const Point&, const TMatrixD&,
                      int = -1, const std::vector<MeasuredHadTauDecayProduct>* = 0);
    MeasuredTauLepton(int, 
                      int, double, double, double, double, 
                      const Point&, double, double,
                      int = -1, const std::vector<MeasuredHadTauDecayProduct>* = 0);
    MeasuredTauLepton(const MeasuredTauLepton&);
    ~MeasuredTauLepton();

    MeasuredTauLepton& 
    operator=(const MeasuredTauLepton&);

    /// return decay type of the tau lepton
    int
    type() const;

    std::string
    type_string() const;

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

    /// return decay mode of the reconstructed hadronic tau decay
    int
    decayMode() const;

    /// return position of tau decay vertex and its uncertainty
    bool
    hasDecayVertex() const;
    const Point&
    decayVertex() const;
    const TMatrixD&
    decayVertexCov() const;
    const TMatrixD&
    decayVertexCovInv() const;
    bool
    decayVertexCovInv_isValid() const;

    /// return individual charged and neutral pions produced in hadronic tau decay
    bool
    hasHadTauDecayProducts() const;
    const std::vector<MeasuredHadTauDecayProduct>&
    hadTauDecayProducts() const;
    const MeasuredHadTauDecayProduct*
    leadChargedHadron() const;

    /// return four-vector in labframe
    const LorentzVector&
    p4() const;

    /// return momentum vector in labframe
    const Vector&
    p3() const;

    /// return flags indicating if tau decayed leptonically or hadronically 
    bool
    isLeptonicTauDecay() const;
    bool
    isHadronicTauDecay() const;

    /// return flag indicating electrons or muons directly originating from LFV Higgs boson decay
    bool
    isPrompt() const;

   protected:
    /// check that MeasuredTauLepton is of valid type
    void
    checkType();

    /// check that mass of MeasuredTauLepton matches given type
    /// and set preciseVisMass data-member
    void
    setMass();

    /// set measuredDecayVertex and covDecayVertex data-members,
    /// compute inverse of covariance matrix
    void
    setDecayVertex(const Point&, const TMatrixD&);
    void
    setDecayVertex(const Point&, double, double);

    /// set measuredHadTauDecayProducts data-member
    void
    setHadTauDecayProducts(const std::vector<MeasuredHadTauDecayProduct>*);

    /// set visible momentum in all coordinates systems
    void
    initialize();

    /// decay type
    int type_;
    std::string type_string_;

    /// charge
    int charge_;

    /// visible momentum in labframe (in polar coordinates)
    double pt_;
    double eta_;
    double phi_;
    double mass_;

    /// visible momentum in labframe (in cartesian coordinates)
    double energy_;
    double px_;
    double py_;
    double pz_;

    /// magnitude of visible momentum in labframe (magnitude)
    double p_;

    /// decay mode (hadronic tau decays only)
    int decayMode_;

    /// position of tau decay vertex and its uncertainty
    bool hasDecayVertex_;
    Point decayVertex_;
    TMatrixD decayVertexCov_;
    TMatrixD decayVertexCovInv_;
    bool decayVertexCovInv_isValid_;

    /// individual charged and neutral pions (hadronic tau decays only)
    bool hasHadTauDecayProducts_;
    std::vector<MeasuredHadTauDecayProduct> hadTauDecayProducts_;
    const MeasuredHadTauDecayProduct* leadChargedHadron_;

    /// visible momentum in labframe (four-vector)
    LorentzVector p4_;

    /// visible momentum in labframe
    Vector p3_;

    /// mass of visible tau decay products (recomputed to reduce rounding errors)
    double preciseVisMass_;

    /// flags indicating if tau decayed leptonically or hadronically
    bool isLeptonicTauDecay_;
    bool isHadronicTauDecay_;

    /// flag indicating electrons or muons directly originating from LFV Higgs boson decay
    bool isPrompt_;
  };

  // auxiliary function for printing MeasuredTauLepton objects (for debugging purposes)
  std::ostream&
  operator<<(std::ostream& os, const std::vector<MeasuredTauLepton>& measuredTauLeptons);

  // auxiliary class for sorting MeasuredTauLepton objects
  struct sortMeasuredTauLeptons
  {
    bool operator() (const MeasuredTauLepton& measuredTauLepton1, const MeasuredTauLepton& measuredTauLepton2);
  };
}

#endif
