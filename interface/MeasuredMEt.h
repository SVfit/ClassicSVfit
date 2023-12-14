#ifndef TauAnalysis_ClassicSVfit_MeasuredMEt_h
#define TauAnalysis_ClassicSVfit_MeasuredMEt_h

#include <TMatrixD.h>

#include <vector>

namespace classic_svFit
{
  class MeasuredMEt
  {
   public:
    /// default constructor
    MeasuredMEt();

    /// constructor for proton-proton collisions
    MeasuredMEt(double, double, const TMatrixD&);

    /// constructor for electron-positron collisions
    MeasuredMEt(double, double, double, double, const TMatrixD&);

    // copy constructor
    MeasuredMEt(const MeasuredMEt&);

    // destructor
    ~MeasuredMEt();

    /// return type of event (either proton-proton or electron-positron collisions)
    int
    type() const;

    /// return missing momentum (Pz and energy components defined for electron-positron collisions only)
    double
    px() const;
    double
    py() const;
    double
    pz() const;
    double
    energy() const;

    /// return uncertainty on missing momentum (covariance matrix and its inverse)
    TMatrixD&
    cov() const;
    TMatrixD&
    covInv() const;
    bool
    covInv_isValid() const;

    /// return constant factor for evaluating ClassicSVfitIntegrand::EvalMEtTF() function
    /// (pre-computed to reduce computing time)
    double
    const_MET() const;

   protected:
    // set cov data-member
    void
    setCov();

    /// type of MeasuredMEt object (either proton-proton or electron-positron collisions)
    int type_;
    std::string type_string_;

    /// missing momentum (pz and energy components defined for electron-positron collisions only)
    double px_;
    double py_;
    double pz_;
    double energy_;

    /// uncertainty on missing momentum (covariance matrix and its inverse)
    TMatrixD cov_;
    TMatrixD covInv_;
    bool covInv_isValid_;

    /// constant factor for evaluating ClassicSVfitIntegrand::EvalMEtTF() function
    /// (pre-computed to reduce computing time)
    double const_MET_;
  };
}

#endif
