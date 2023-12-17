#ifndef TauAnalysis_ClassicSVfit_MeasuredEvent_h
#define TauAnalysis_ClassicSVfit_MeasuredEvent_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredMEt.h"       // MeasuredMEt
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h" // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // kEventType

#include <TMatrixD.h>                                             // TMatrixD

#include <vector>                                                 // std::vector<>

namespace classic_svFit
{
  class MeasuredEvent
  {
   public:
    MeasuredEvent();
    MeasuredEvent(const std::vector<MeasuredTauLepton>&, const std::vector<MeasuredMEt>&);
    MeasuredEvent(const std::vector<MeasuredTauLepton>&, const std::vector<MeasuredMEt>&,
                  const Point&, const TMatrixD&);
    MeasuredEvent(const MeasuredEvent&);
    ~MeasuredEvent();

    MeasuredEvent& 
    operator=(const MeasuredEvent&);

    /// return type of event (either proton-proton or electron-positron collisions)
    int
    type() const;

    /// return visible tau decay products
    const std::vector<MeasuredTauLepton>&
    tauLeptons() const;

    const MeasuredTauLepton&
    tauPlus() const;
    const MeasuredTauLepton&
    tauMinus() const;

    /// return missing momentum
    const std::vector<MeasuredMEt>&
    MEt() const;

    /// return position of primary vertex and its uncertainty
    bool
    hasPrimaryVertex() const;
    const Point&
    primaryVertex() const;
    const TMatrixD&
    primaryVertexCov() const;
    const TMatrixD&
    primaryVertexCovInv() const;
    bool
    primaryVertexCovInv_isValid() const;

   protected:
    /// set measuredTauLeptons data-member
    void
    setTauLeptons(const std::vector<MeasuredTauLepton>&);

    /// set measuredMEt data-member
    void
    setMEt(const std::vector<MeasuredMEt>&);

    /// set measuredPrimaryVertex and covPrimaryVertex data-members,
    /// compute inverse of covariance matrix
    void
    setPrimaryVertex(const Point&, const TMatrixD&);

    /// type of event (either proton-proton or electron-positron collisions)
    int type_;

    /// visible tau decay products
    std::vector<MeasuredTauLepton> tauLeptons_;
    const MeasuredTauLepton* tauPlus_;
    const MeasuredTauLepton* tauMinus_;

    /// missing momentum
    std::vector<MeasuredMEt> MEt_;

    /// position of tau primary vertex and its uncertainty
    bool hasPrimaryVertex_;
    Point primaryVertex_;
    TMatrixD primaryVertexCov_;
    TMatrixD primaryVertexCovInv_;
    bool primaryVertexCovInv_isValid_;
  };
}

#endif
