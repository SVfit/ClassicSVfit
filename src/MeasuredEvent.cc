#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"

#include <algorithm> // std::sort()
#include <assert.h>  // assert()
#include <iostream>  // std::cerr, std::cout, std::endl

namespace classic_svFit
{

MeasuredEvent::MeasuredEvent()
  : type_(kUndefinedCollisionType)
  , tauPlus_(nullptr)
  , tauMinus_(nullptr)
  , hasPrimaryVertex_(false)
  , primaryVertexCovInv_isValid_(false)
{}

MeasuredEvent::MeasuredEvent(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const std::vector<MeasuredMEt>& measuredMEt)
  : tauPlus_(nullptr)
  , tauMinus_(nullptr)
  , hasPrimaryVertex_(false)
  , primaryVertexCovInv_isValid_(false)
{
  setTauLeptons(measuredTauLeptons);
  setMEt(measuredMEt);
}

MeasuredEvent::MeasuredEvent(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const std::vector<MeasuredMEt>& measuredMEt,
                             const Point& primaryVertex, const TMatrixD& primaryVertexCov)
  : tauPlus_(nullptr)
  , tauMinus_(nullptr)
  , hasPrimaryVertex_(true)
  , primaryVertexCovInv_isValid_(false)
{
  setTauLeptons(measuredTauLeptons);
  setMEt(measuredMEt);
  setPrimaryVertex(primaryVertex, primaryVertexCov);
}

MeasuredEvent::MeasuredEvent(const MeasuredEvent& measuredEvent)
  : type_(measuredEvent.type_)
  , tauLeptons_(measuredEvent.tauLeptons_)
  , tauPlus_(measuredEvent.tauPlus_)
  , tauMinus_(measuredEvent.tauMinus_)
  , MEt_(measuredEvent.MEt_)
  , hasPrimaryVertex_(measuredEvent.hasPrimaryVertex_)
  , primaryVertex_(measuredEvent.primaryVertex_)
  , primaryVertexCov_(measuredEvent.primaryVertexCov_)
  , primaryVertexCovInv_(measuredEvent.primaryVertexCovInv_)
  , primaryVertexCovInv_isValid_(measuredEvent.primaryVertexCovInv_isValid_)
{}

MeasuredEvent::~MeasuredEvent()
{}

MeasuredEvent& 
MeasuredEvent::operator=(const MeasuredEvent& measuredEvent)
{
  type_ = measuredEvent.type_;
  tauLeptons_ = measuredEvent.tauLeptons_;
  tauPlus_ = measuredEvent.tauPlus_;
  tauMinus_ = measuredEvent.tauMinus_;
  MEt_ = measuredEvent.MEt_;
  hasPrimaryVertex_ = measuredEvent.hasPrimaryVertex_;
  primaryVertex_ = measuredEvent.primaryVertex_;
  primaryVertexCov_.ResizeTo(measuredEvent.primaryVertexCov_.GetNrows(),measuredEvent.primaryVertexCov_.GetNcols());
  primaryVertexCov_ = measuredEvent.primaryVertexCov_;
  primaryVertexCovInv_.ResizeTo(measuredEvent.primaryVertexCovInv_.GetNrows(),measuredEvent.primaryVertexCovInv_.GetNcols());
  primaryVertexCovInv_ = measuredEvent.primaryVertexCovInv_;
  primaryVertexCovInv_isValid_ = measuredEvent.primaryVertexCovInv_isValid_;
  return *this;
}

int
MeasuredEvent::type() const
{
  return type_;
}

const std::vector<MeasuredTauLepton>&
MeasuredEvent::tauLeptons() const
{
  return tauLeptons_;
}

const MeasuredTauLepton&
MeasuredEvent::tauPlus() const
{
  assert(tauPlus_);
  return *tauPlus_;
}

const MeasuredTauLepton&
MeasuredEvent::tauMinus() const
{
  assert(tauMinus_);
  return *tauMinus_;
}

const std::vector<MeasuredMEt>&
MeasuredEvent::MEt() const
{
  return MEt_;
}

bool
MeasuredEvent::hasPrimaryVertex() const
{
  return hasPrimaryVertex_;
}

const Point&
MeasuredEvent::primaryVertex() const
{
  return primaryVertex_;
}

const TMatrixD&
MeasuredEvent::primaryVertexCov() const
{
  return primaryVertexCov_;
}

const TMatrixD&
MeasuredEvent::primaryVertexCovInv() const
{
  return primaryVertexCovInv_;
}

bool
MeasuredEvent::primaryVertexCovInv_isValid() const
{
  return primaryVertexCovInv_isValid_;
}

void
MeasuredEvent::setTauLeptons(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  tauLeptons_ = measuredTauLeptons;
  std::sort(tauLeptons_.begin(), tauLeptons_.end(), sortMeasuredTauLeptons());
  for ( const MeasuredTauLepton& tauLepton : tauLeptons_ )
  {
    if      ( tauLepton.charge() > 0 ) tauPlus_  = &tauLepton;
    else if ( tauLepton.charge() < 0 ) tauMinus_ = &tauLepton;
  }
  assert(tauPlus_ && tauMinus_);
}

void
MeasuredEvent::setMEt(const std::vector<MeasuredMEt>& measuredMEt)
{
  if ( measuredMEt.size() == 0 )
  {
    std::cerr << "ERROR: No MET given !!" << std::endl;
    assert(0);
  }
  type_ = measuredMEt.front().type();
  for ( const MeasuredMEt& measuredMEt_i : measuredMEt )
  {
    if ( measuredMEt_i.type() != type_ )
    {
      std::cerr << "ERROR: Not all MET objects are of the same type !!" << std::endl;
      assert(0);
    }
  }
  MEt_ = measuredMEt;
}

void
MeasuredEvent::setPrimaryVertex(const Point& primaryVertex, const TMatrixD& primaryVertexCov)
{
  double primaryVertexX = roundToNdigits(primaryVertex.x());
  double primaryVertexY = roundToNdigits(primaryVertex.y());
  double primaryVertexZ = roundToNdigits(primaryVertex.z());
  primaryVertex_ = Point(primaryVertexX, primaryVertexY, primaryVertexZ);

  assert(primaryVertexCov.GetNrows() == primaryVertexCov.GetNcols());
  int dim = primaryVertexCov.GetNrows();
  primaryVertexCov_.ResizeTo(dim,dim);
  for ( int iRow = 0; iRow < dim; ++iRow )
  {
    for ( int iColumn = 0; iColumn < dim; ++iColumn )
    {
      primaryVertexCov_(iRow,iColumn) = roundToNdigits(primaryVertexCov(iRow,iColumn));
    }
  }
  //std::cout << "primaryVertexCov:" << std::endl;
  //primaryVertexCov_.Print();
  if ( primaryVertexCov_.Determinant() == 0. )
  {
    std::cout << "primaryVertexCov:" << std::endl;
    primaryVertexCov_.Print();
    std::cerr << "ERROR: Failed to invert matrix primaryVertexCov (det=0) !!" << std::endl;
    return;
  }

  primaryVertexCovInv_.ResizeTo(dim,dim);
  primaryVertexCovInv_ = TMatrixD(TMatrixD::kInverted, primaryVertexCov_);
  primaryVertexCovInv_isValid_ = true;
  //std::cout << "primaryVertexCovInv:" << std::endl;
  //primaryVertexCovInv_.Print();
}

}
