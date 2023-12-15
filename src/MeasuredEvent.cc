#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"

#include <algorithm> // std::sort()
#include <assert.h>  // assert()
#include <iostream>  // std::cerr, std::cout, std::endl

namespace classic_svFit
{

MeasuredEvent::MeasuredEvent()
  : type_(kUndefinedCollisionType)
  , measuredTauPlus_(nullptr)
  , measuredTauMinus_(nullptr)
  , hasPrimaryVertex_(false)
  , covInvPrimaryVertex_isValid_(false)
{}

MeasuredEvent::MeasuredEvent(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const std::vector<MeasuredMEt>& measuredMEt)
  : measuredTauPlus_(nullptr)
  , measuredTauMinus_(nullptr)
  , hasPrimaryVertex_(false)
  , covInvPrimaryVertex_isValid_(false)
{
  setTauLeptons(measuredTauLeptons);
  setMEt(measuredMEt);
}

MeasuredEvent::MeasuredEvent(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const std::vector<MeasuredMEt>& measuredMEt,
                             const Point& measuredPrimaryVertex, const TMatrixD& covPrimaryVertex)
  : measuredTauPlus_(nullptr)
  , measuredTauMinus_(nullptr)
  , hasPrimaryVertex_(true)
  , covInvPrimaryVertex_isValid_(false)
{
  setTauLeptons(measuredTauLeptons);
  setMEt(measuredMEt);
  setPrimaryVertex(measuredPrimaryVertex, covPrimaryVertex);
}

MeasuredEvent::MeasuredEvent(const MeasuredEvent& measuredEvent)
  : type_(measuredEvent.type_)
  , measuredTauLeptons_(measuredEvent.measuredTauLeptons_)
  , measuredTauPlus_(measuredEvent.measuredTauPlus_)
  , measuredTauMinus_(measuredEvent.measuredTauMinus_)
  , measuredMEt_(measuredEvent.measuredMEt_)
  , hasPrimaryVertex_(measuredEvent.hasPrimaryVertex_)
  , measuredPrimaryVertex_(measuredEvent.measuredPrimaryVertex_)
  , covPrimaryVertex_(measuredEvent.covPrimaryVertex_)
  , covInvPrimaryVertex_(measuredEvent.covInvPrimaryVertex_)
  , covInvPrimaryVertex_isValid_(measuredEvent.covInvPrimaryVertex_isValid_)
{}

MeasuredEvent::~MeasuredEvent()
{}

MeasuredEvent& 
MeasuredEvent::operator=(const MeasuredEvent& measuredEvent)
{
  type_ = measuredEvent.type_;
  measuredTauLeptons_ = measuredEvent.measuredTauLeptons_;
  measuredTauPlus_ = measuredEvent.measuredTauPlus_;
  measuredTauMinus_ = measuredEvent.measuredTauMinus_;
  measuredMEt_ = measuredEvent.measuredMEt_;
  hasPrimaryVertex_ = measuredEvent.hasPrimaryVertex_;
  measuredPrimaryVertex_ = measuredEvent.measuredPrimaryVertex_;
  covPrimaryVertex_.ResizeTo(measuredEvent.covPrimaryVertex_.GetNrows(),measuredEvent.covPrimaryVertex_.GetNcols());
  covPrimaryVertex_ = measuredEvent.covPrimaryVertex_;
  covInvPrimaryVertex_.ResizeTo(measuredEvent.covInvPrimaryVertex_.GetNrows(),measuredEvent.covInvPrimaryVertex_.GetNcols());
  covInvPrimaryVertex_ = measuredEvent.covInvPrimaryVertex_;
  covInvPrimaryVertex_isValid_ = measuredEvent.covInvPrimaryVertex_isValid_;
  return *this;
}

int
MeasuredEvent::type() const
{
  return type_;
}

const std::vector<MeasuredTauLepton>&
MeasuredEvent::measuredTauLeptons() const
{
  return measuredTauLeptons_;
}

const MeasuredTauLepton&
MeasuredEvent::measuredTauPlus() const
{
  assert(measuredTauPlus_);
  return *measuredTauPlus_;
}

const MeasuredTauLepton&
MeasuredEvent::measuredTauMinus() const
{
  assert(measuredTauMinus_);
  return *measuredTauMinus_;
}

const std::vector<MeasuredMEt>&
MeasuredEvent::measuredMEt() const
{
  return measuredMEt_;
}

bool
MeasuredEvent::hasPrimaryVertex() const
{
  return hasPrimaryVertex_;
}

const Point&
MeasuredEvent::measuredPrimaryVertex() const
{
  return measuredPrimaryVertex_;
}

const TMatrixD&
MeasuredEvent::covPrimaryVertex() const
{
  return covPrimaryVertex_;
}

const TMatrixD&
MeasuredEvent::covInvPrimaryVertex() const
{
  return covInvPrimaryVertex_;
}

bool
MeasuredEvent::covInvPrimaryVertex_isValid() const
{
  return covInvPrimaryVertex_isValid_;
}

void
MeasuredEvent::setTauLeptons(const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  measuredTauLeptons_ = measuredTauLeptons;
  std::sort(measuredTauLeptons_.begin(), measuredTauLeptons_.end(), sortMeasuredTauLeptons());
  for ( const MeasuredTauLepton& measuredTauLepton : measuredTauLeptons_ )
  {
    if      ( measuredTauLepton.charge() > 0 ) measuredTauPlus_ = &measuredTauLepton;
    else if ( measuredTauLepton.charge() < 0 ) measuredTauMinus_ = &measuredTauLepton;
  }
  assert(measuredTauPlus_ && measuredTauMinus_);
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
  measuredMEt_ = measuredMEt;
}

void
MeasuredEvent::setPrimaryVertex(const Point& measuredPrimaryVertex, const TMatrixD& covPrimaryVertex)
{
  double measuredPrimaryVertexX = roundToNdigits(measuredPrimaryVertex.x());
  double measuredPrimaryVertexY = roundToNdigits(measuredPrimaryVertex.y());
  double measuredPrimaryVertexZ = roundToNdigits(measuredPrimaryVertex.z());
  measuredPrimaryVertex_ = Point(measuredPrimaryVertexX, measuredPrimaryVertexY, measuredPrimaryVertexZ);

  assert(covPrimaryVertex.GetNrows() == covPrimaryVertex.GetNcols());
  int dim = covPrimaryVertex.GetNrows();
  covPrimaryVertex_.ResizeTo(dim,dim);
  for ( int iRow = 0; iRow < dim; ++iRow )
  {
    for ( int iColumn = 0; iColumn < dim; ++iColumn )
    {
      covPrimaryVertex_(iRow,iColumn) = roundToNdigits(covPrimaryVertex(iRow,iColumn));
    }
  }
  if ( covPrimaryVertex_.Determinant() == 0. )
  {
    std::cout << "covPrimaryVertex:" << std::endl;
    covPrimaryVertex_.Print();
    std::cerr << "ERROR: Failed to invert matrix covPrimaryVertex (det=0) !!" << std::endl;
    return;
  }
  covInvPrimaryVertex_.ResizeTo(dim,dim);
  covInvPrimaryVertex_ = TMatrixD(TMatrixD::kInverted, covPrimaryVertex_);
  covInvPrimaryVertex_isValid_ = true;
}

}
