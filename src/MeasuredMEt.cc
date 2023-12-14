#include "TauAnalysis/ClassicSVfit/interface/MeasuredMEt.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // kEventType

#include <TMath.h>                                                // TMath::Pi()
#include <TMatrixDSym.h>                                          // TMatrixDSym
#include <TMatrixDSymEigen.h>                                     // TMatrixDSymEigen
#include <TVectorD.h>                                             // TVectorD

#include <assert.h>                                               // assert()
#include <cmath>                                                  // atan2(), pow(), std::sqrt()
#include <iostream>                                               // std::cout, std::endl

namespace classic_svFit
{

MeasuredMEt::MeasuredMEt()
  : type_(kUndefinedCollisionType)
  , type_string_("undefined")
  , px_(0.)
  , py_(0.)
  , pz_(0.)
  , energy_(0.)
  , covInv_isValid_(false)
  , const_MET_(0.)
{}

MeasuredMEt::MeasuredMEt(double px, double py, const TMatrixD& cov)
  : type_(kProtonProtonCollisions)
  , type_string_("proton-proton")
  , px_(px)
  , py_(py)
  , pz_(0.)
  , energy_(0.)
  , covInv_isValid_(false)
  , const_MET_(0.)
{
  setCov(cov);
}

MeasuredMEt::MeasuredMEt(double px, double py, double pz, double energy, const TMatrixD& cov)
  : type_(kElectronPositronCollisions)
  , type_string_("electron-positron")
  , px_(px)
  , py_(py)
  , pz_(pz)
  , energy_(energy)
  , covInv_isValid_(false)
  , const_MET_(0.)
{
  setCov(cov);
}

MeasuredMEt::MeasuredMEt(const MeasuredMEt& measuredMEt)
  : type_(measuredMEt.type_)
  , type_string_(measuredMEt.type_string_)
  , px_(measuredMEt.px_)
  , py_(measuredMEt.py_)
  , pz_(measuredMEt.pz_)
  , energy_(measuredMEt.energy_)
  , cov_(measuredMEt.cov_)
  , covInv_(measuredMEt.covInv_)
  , covInv_isValid_(measuredMEt.covInv_isValid_)
  , const_MET_(measuredMEt.const_MET_)
{}

MeasuredMEt::~MeasuredMEt()
{}

int
MeasuredMEt::type() const
{
  return type_;
}

double
MeasuredMEt::px() const
{
  return px_;
}

double
MeasuredMEt::py() const
{
  return py_;
}
 
double
MeasuredMEt::pz() const
{
  return pz_;
}

double
MeasuredMEt::energy() const
{
  return energy_;
}

const TMatrixD&
MeasuredMEt::cov() const
{
  return cov_;
}

const TMatrixD&
MeasuredMEt::covInv() const
{
  return covInv_;
}

bool
MeasuredMEt::covInv_isValid() const
{
  return covInv_isValid_;
}

double
MeasuredMEt::const_MET() const
{
  return const_MET_;
}

void
MeasuredMEt::setCov(const TMatrixD& cov)
{
  assert(cov.GetNrows() == cov.GetNcols());
  int dim = cov.GetNrows();
  if ( (type_ == kProtonProtonCollisions     && dim != 2) ||
       (type_ == kElectronPositronCollisions && dim != 4) )
  {
    std::cerr << "ERROR: Size of covariance matrix = (" << cov.GetNrows() << "," << cov.GetNcols() << ")" 
              << " does not match size expected for " << type_string_ << " collisions !!" << std::endl;
    assert(0);
  }
  cov_.ResizeTo(dim,dim);
  for ( int iRow = 0; iRow < dim; ++iRow )
  {
    for ( int iColumn = 0; iColumn < dim; ++iColumn )
    {
      cov_(iRow,iColumn) = roundToNdigits(cov(iRow,iColumn));
    }
  }
  if ( cov_.Determinant() == 0. )
  {
    std::cout << "cov:" << std::endl;
    cov_.Print();
    std::cerr << "ERROR: Failed to invert matrix cov (det=0) !!" << std::endl;
    return;
  }
  covInv_.ResizeTo(dim,dim);
  covInv_ = TMatrixD(TMatrixD::kInverted, cov_);
  covInv_isValid_ = true;
  const_MET_ = 1./(pow(2.*TMath::Pi(), 0.5*dim)*std::sqrt(cov_.Determinant()));
}

//---------------------------------------------------------------------------------------------------
// auxiliary function for printing MeasuredMEt object (for debugging purposes)
std::ostream&
operator<<(std::ostream& os, const MeasuredMEt& measuredMEt)
{
  os << "MET: Px = " << measuredMEt.px() << ", Py = " << measuredMEt.py();
  if ( measuredMEt.type() == kElectronPositronCollisions )
  {
    os << ", Pz = " << measuredMEt.pz() << ", E = " << measuredMEt.energy();
  }
  os << std::endl;
  const TMatrixD& cov = measuredMEt.cov();
  os << "cov:" << std::endl;
  cov.Print();
  if ( measuredMEt.type() == kProtonProtonCollisions )
  {
    assert(cov.GetNrows() == 2 && cov.GetNcols() == 2);
    int dim = cov.GetNrows();
    TMatrixDSym cov_sym(dim);
    for ( int iRow = 0; iRow < dim; ++iRow )
    {
      for ( int iColumn = 0; iColumn < dim; ++iColumn )
      {
        cov_sym(iRow,iColumn) = cov(iRow,iColumn);
      }
    }
    TMatrixD EigenVectors(dim,dim);
    EigenVectors = TMatrixDSymEigen(cov_sym).GetEigenVectors();
    os << "Eigenvectors =  { ";
    for ( int iEigenVector = 0; iEigenVector < dim; ++iEigenVector )
    {
      if ( iEigenVector > 0 ) std::cout << ", ";
      os << "{ " << EigenVectors(0,iEigenVector) << ", " << EigenVectors(1,iEigenVector) 
         << " (phi = " << atan2(EigenVectors(1,iEigenVector), EigenVectors(0,iEigenVector)) << ") }";
    }
    os << std::endl;
    TVectorD EigenValues(dim);
    EigenValues = TMatrixDSymEigen(cov_sym).GetEigenValues();
    os << "sqrt(Eigenvalues) = ";
    for ( int iElement = 0; iElement < dim; ++iElement )
    {
      if ( iElement > 0 ) std::cout << ", ";
      os << std::sqrt(EigenValues(iElement));
    }
    os << std::endl;
  }
  return os;
}
//---------------------------------------------------------------------------------------------------

}
