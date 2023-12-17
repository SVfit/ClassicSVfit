#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

#include <algorithm> // std::sort()
#include <assert.h>  // assert()
#include <cmath>     // cos(), cosh(), sin(), sinh(), std::sqrt()
#include <iostream>  // std::cerr, std::cout, std::endl

namespace classic_svFit
{

MeasuredTauLepton::MeasuredTauLepton()
  : type_(kUndefinedDecayType)
  , charge_(0)
  , pt_(0.)
  , eta_(0.)
  , phi_(0.)
  , mass_(0.)
  , decayMode_(-1)
  , hasDecayVertex_(false)
  , decayVertexCovInv_isValid_(false)
  , hasHadTauDecayProducts_(false)
{
  initialize();
}

MeasuredTauLepton::MeasuredTauLepton(int type, int charge, double pt, double eta, double phi, double mass,
                                     int decayMode, const std::vector<MeasuredHadTauDecayProduct>* hadTauDecayProducts)
  : type_(type)
  , type_string_("undefined")
  , charge_(charge)
  , pt_(roundToNdigits(pt))
  , eta_(roundToNdigits(eta))
  , phi_(roundToNdigits(phi))
  , mass_(roundToNdigits(mass))
  , decayMode_(decayMode)
  , hasDecayVertex_(false)
  , decayVertexCovInv_isValid_(false)
  , hasHadTauDecayProducts_(hadTauDecayProducts != nullptr) 
{
  checkType();
  setMass();
  setHadTauDecayProducts(hadTauDecayProducts);
  initialize();
}

MeasuredTauLepton::MeasuredTauLepton(int type, 
                                     int charge, double pt, double eta, double phi, double mass, 
                                     const Point& decayVertex, const TMatrixD& decayVertexCov,
                                     int decayMode, const std::vector<MeasuredHadTauDecayProduct>* hadTauDecayProducts)
  : type_(type)
  , type_string_("undefined")
  , charge_(charge)
  , pt_(roundToNdigits(pt))
  , eta_(roundToNdigits(eta))
  , phi_(roundToNdigits(phi))
  , mass_(roundToNdigits(mass))
  , decayMode_(decayMode)
  , hasDecayVertex_(true)
  , decayVertexCovInv_isValid_(false)
  , hasHadTauDecayProducts_(hadTauDecayProducts != nullptr)
  , leadChargedHadron_(nullptr)
{
  checkType();
  setMass();
  setHadTauDecayProducts(hadTauDecayProducts);
  setDecayVertex(decayVertex, decayVertexCov);
  initialize();
}

MeasuredTauLepton::MeasuredTauLepton(int type, 
                                     int charge, double pt, double eta, double phi, double mass, 
                                     const Point& decayVertex, double decayVertexSigma_parl, double decayVertexSigma_perp,
                                     int decayMode, const std::vector<MeasuredHadTauDecayProduct>* hadTauDecayProducts)
  : type_(type)
  , type_string_("undefined")
  , charge_(charge)
  , pt_(roundToNdigits(pt))
  , eta_(roundToNdigits(eta))
  , phi_(roundToNdigits(phi))
  , mass_(roundToNdigits(mass))
  , decayMode_(decayMode)
  , hasDecayVertex_(true)
  , decayVertexCovInv_isValid_(false)
  , hasHadTauDecayProducts_(hadTauDecayProducts != nullptr)
  , leadChargedHadron_(nullptr)
{
  checkType();
  setMass();
  setHadTauDecayProducts(hadTauDecayProducts);
  setDecayVertex(decayVertex, decayVertexSigma_parl, decayVertexSigma_perp);
  initialize();
}

MeasuredTauLepton::MeasuredTauLepton(const MeasuredTauLepton& measuredTauLepton)
  : type_(measuredTauLepton.type_)
  , charge_(measuredTauLepton.charge_)
  , pt_(measuredTauLepton.pt_)
  , eta_(measuredTauLepton.eta_)
  , phi_(measuredTauLepton.phi_)
  , mass_(measuredTauLepton.mass_)
  , decayMode_(measuredTauLepton.decayMode_)
  , hasDecayVertex_(measuredTauLepton.hasDecayVertex_)
  , decayVertex_(measuredTauLepton.decayVertex_)
  , decayVertexCov_(measuredTauLepton.decayVertexCov_)
  , decayVertexCovInv_(measuredTauLepton.decayVertexCovInv_)
  , decayVertexCovInv_isValid_(measuredTauLepton.decayVertexCovInv_isValid_)
  , hasHadTauDecayProducts_(measuredTauLepton.hasHadTauDecayProducts_)
{
  preciseVisMass_ = measuredTauLepton.mass();
  setHadTauDecayProducts(&measuredTauLepton.hadTauDecayProducts_);
  initialize();
}

MeasuredTauLepton::~MeasuredTauLepton()
{}

MeasuredTauLepton& 
MeasuredTauLepton::operator=(const MeasuredTauLepton& measuredTauLepton)
{
  type_ = measuredTauLepton.type_;
  charge_ = measuredTauLepton.charge_;
  pt_ = measuredTauLepton.pt_;
  eta_ = measuredTauLepton.eta_;
  phi_ = measuredTauLepton.phi_;
  mass_ = measuredTauLepton.mass_;
  decayMode_ = measuredTauLepton.decayMode_;
  hasDecayVertex_ = measuredTauLepton.hasDecayVertex_;
  decayVertex_ = measuredTauLepton.decayVertex_;
  decayVertexCov_.ResizeTo(measuredTauLepton.decayVertexCov_.GetNrows(),measuredTauLepton.decayVertexCov_.GetNcols());
  decayVertexCov_ = measuredTauLepton.decayVertexCov_;
  decayVertexCovInv_.ResizeTo(measuredTauLepton.decayVertexCovInv_.GetNrows(),measuredTauLepton.decayVertexCovInv_.GetNcols());
  decayVertexCovInv_ = measuredTauLepton.decayVertexCovInv_;
  decayVertexCovInv_isValid_ = measuredTauLepton.decayVertexCovInv_isValid_;
  hasHadTauDecayProducts_ = measuredTauLepton.hasHadTauDecayProducts_;
  preciseVisMass_ = measuredTauLepton.mass();
  setHadTauDecayProducts(&measuredTauLepton.hadTauDecayProducts_);
  initialize();
  return *this;
}

int
MeasuredTauLepton::type() const 
{ 
  return type_; 
}

std::string
MeasuredTauLepton::type_string() const
{
  return type_string_; 
}

int
MeasuredTauLepton::charge() const 
{ 
  return charge_; 
}

double
MeasuredTauLepton::pt() const 
{ 
  return pt_; 
}

double
MeasuredTauLepton::eta() const 
{ 
  return eta_; 
}

double
MeasuredTauLepton::phi() const 
{ 
  return phi_; 
}

double
MeasuredTauLepton::mass() const 
{ 
  return preciseVisMass_; 
}

double
MeasuredTauLepton::energy() const 
{ 
  return energy_; 
}

double
MeasuredTauLepton::px() const 
{ 
  return px_; 
}

double
MeasuredTauLepton::py() const 
{ 
  return py_; 
}

double
MeasuredTauLepton::pz() const 
{ 
  return pz_; 
}

double
MeasuredTauLepton::p() const 
{
  return p_; 
}

int
MeasuredTauLepton::decayMode() const 
{
  return decayMode_;
}

bool
MeasuredTauLepton::hasDecayVertex() const
{
  return hasDecayVertex_;
}

const Point&
MeasuredTauLepton::decayVertex() const
{
  return decayVertex_;
}

const TMatrixD&
MeasuredTauLepton::decayVertexCov() const
{
  return decayVertexCov_;
}

const TMatrixD&
MeasuredTauLepton::decayVertexCovInv() const
{
  return decayVertexCovInv_;
}

bool
MeasuredTauLepton::decayVertexCovInv_isValid() const
{
  return decayVertexCovInv_isValid_;
}

bool
MeasuredTauLepton::hasHadTauDecayProducts() const
{
  return hasHadTauDecayProducts_;
}

const std::vector<MeasuredHadTauDecayProduct>&
MeasuredTauLepton::hadTauDecayProducts() const
{
  return hadTauDecayProducts_;
}

const MeasuredHadTauDecayProduct*
MeasuredTauLepton::leadChargedHadron() const
{
  return leadChargedHadron_;
}

const LorentzVector&
MeasuredTauLepton::p4() const 
{
  return p4_; 
}

const Vector&
MeasuredTauLepton::p3() const 
{ 
  return p3_; 
}

void
MeasuredTauLepton::checkType()
{
  if      ( type_ == kPrompt         ) type_string_ = "prompt lepton (lepton-flavor-violating)";
  else if ( type_ == kTauToElecDecay ) type_string_ = "tau -> electron decay";
  else if ( type_ == kTauToMuDecay   ) type_string_ = "tau -> muon decay";
  else if ( type_ == kTauToHadDecay  ) type_string_ = "tau -> had decay";
  else
  {
    std::cerr << "ERROR: Invalid type " << type_ << " declared for leg:" 
              << " Pt = " << pt_ << ", eta = " << eta_ << ", phi = " << phi_ << ", mass = " << mass_ << " !!" << std::endl;
    assert(0);
  }
}

void
MeasuredTauLepton::setMass()
{
  double minVisMass = electronMass;
  double maxVisMass = tauLeptonMass;
  if ( type_ == kTauToElecDecay )
  {
    minVisMass = electronMass;
    maxVisMass = minVisMass;
  } 
  else if ( type_ == kTauToMuDecay ) 
  {
    minVisMass = muonMass;
    maxVisMass = minVisMass;
  }
  else if ( type_ == kTauToHadDecay )
  {
    if ( decayMode_ == -1 )
    {
      minVisMass = chargedPionMass;
      maxVisMass = 1.5;
    } 
    else if ( decayMode_ == 0 )
    {
      minVisMass = chargedPionMass;
      maxVisMass = minVisMass;
    }
    else
    {
      minVisMass = 0.3;
      maxVisMass = 1.5;
    }
  }
  preciseVisMass_ = mass_;
  if ( preciseVisMass_ < (0.9*minVisMass) || preciseVisMass_ > (1.1*maxVisMass) ) 
  {
    std::cerr << "WARNING: " << type_string_;
    if ( type_ == kTauToHadDecay ) std::cerr << " (decayMode = " << decayMode_ << ")";
    std::cerr << " declared for leg:" 
              << " Pt = " << pt_ << ", eta = " << eta_ << ", phi = " << phi_ << ", mass = " << mass_ << " !!" << std::endl;
    std::cerr << " (mass expected within the range = " << minVisMass << ".." << maxVisMass << ")" << std::endl;
  }
  if ( preciseVisMass_ < minVisMass ) preciseVisMass_ = minVisMass;
  if ( preciseVisMass_ > maxVisMass ) preciseVisMass_ = maxVisMass;
}

void
MeasuredTauLepton::setDecayVertex(const Point& decayVertex, const TMatrixD& decayVertexCov)
{
  double decayVertexX = roundToNdigits(decayVertex.x());
  double decayVertexY = roundToNdigits(decayVertex.y());
  double decayVertexZ = roundToNdigits(decayVertex.z());
  decayVertex_ = Point(decayVertexX, decayVertexY, decayVertexZ);

  assert(decayVertexCov.GetNrows() == decayVertexCov.GetNcols());
  int dim = decayVertexCov.GetNrows();
  decayVertexCov_.ResizeTo(dim,dim);
  for ( int iRow = 0; iRow < dim; ++iRow )
  {
    for ( int iColumn = 0; iColumn < dim; ++iColumn )
    {
      decayVertexCov_(iRow,iColumn) = roundToNdigits(decayVertexCov(iRow,iColumn));
    }
  }
  //std::cout << "decayVertexCov:" << std::endl;
  //decayVertexCov_.Print();

  if ( !leadChargedHadron_ )
  {
    std::cerr << "ERROR: Failed to find leading charged hadron !!" << std::endl;
    assert(0);
  }

  Vector r, n, k;
  get_localCoordinateSystem(leadChargedHadron_->p3(), r, n, k);

  TMatrixD rotMatrix_xyz2rnk = get_rotationMatrix(r, n, k);
  TMatrixD rotMatrix_rnk2xyz = get_rotationMatrixInv(r, n, k);

  TMatrixD decayVertexCov_local = rotateCovMatrix(decayVertexCov_, rotMatrix_xyz2rnk);
  //std::cout << "decayVertexCov(local):" << std::endl;
  //decayVertexCov_local.Print();
  if ( decayVertexCov_local.Determinant() == 0. )
  {
    std::cout << "decayVertexCov_local:" << std::endl;
    decayVertexCov_local.Print();
    std::cerr << "ERROR: Failed to invert matrix decayVertexCov_local (det=0) !!" << std::endl;
    return;
  }

  TMatrixD decayVertexCovInv_local(TMatrixD::kInverted, decayVertexCov_local);
  decayVertexCovInv_.ResizeTo(dim,dim);
  decayVertexCovInv_ = rotateCovMatrix(decayVertexCovInv_local, rotMatrix_rnk2xyz);
  decayVertexCovInv_isValid_ = true;
}

void
MeasuredTauLepton::setDecayVertex(const Point& decayVertex, double decayVertexSigma_parl, double decayVertexSigma_perp)
{
  double decayVertexX = roundToNdigits(decayVertex.x());
  double decayVertexY = roundToNdigits(decayVertex.y());
  double decayVertexZ = roundToNdigits(decayVertex.z());
  decayVertex_ = Point(decayVertexX, decayVertexY, decayVertexZ);

  if ( !leadChargedHadron_ )
  {
    std::cerr << "ERROR: Failed to find leading charged hadron !!" << std::endl;
    assert(0);
  }

  Vector r, n, k;
  get_localCoordinateSystem(leadChargedHadron_->p3(), r, n, k);

  TMatrixD decayVertexCov_local(3,3);
  decayVertexCov_local(0,0) = square(decayVertexSigma_perp);
  decayVertexCov_local(1,1) = square(decayVertexSigma_perp);
  decayVertexCov_local(2,2) = square(decayVertexSigma_parl);
  //std::cout << "decayVertexCov(local):" << std::endl;
  //decayVertexCov_local.Print();

  TMatrixD decayVertexCovInv_local(3,3);
  decayVertexCovInv_local(0,0) = 1./square(decayVertexSigma_perp);
  decayVertexCovInv_local(1,1) = 1./square(decayVertexSigma_perp);
  decayVertexCovInv_local(2,2) = 1./square(decayVertexSigma_parl);

  TMatrixD rotMatrix_rnk2xyz = get_rotationMatrixInv(r, n, k);

  decayVertexCov_.ResizeTo(3,3);
  decayVertexCov_ = rotateCovMatrix(decayVertexCov_local, rotMatrix_rnk2xyz);
  //std::cout << "decayVertexCov:" << std::endl;
  //decayVertexCov_.Print();

  decayVertexCovInv_.ResizeTo(3,3);
  decayVertexCovInv_ = rotateCovMatrix(decayVertexCovInv_local, rotMatrix_rnk2xyz);
  decayVertexCovInv_isValid_ = true;
}

void
MeasuredTauLepton::setHadTauDecayProducts(const std::vector<MeasuredHadTauDecayProduct>* hadTauDecayProducts)
{
  if ( hadTauDecayProducts )
  {
    if ( type_ != kTauToHadDecay )
    {
      std::cerr << "ERROR: " << type_string_ << " declared for leg, but hadronic tau decay products given !!" << std::endl;
      assert(0);
    }
    
    hadTauDecayProducts_ = (*hadTauDecayProducts);
    std::sort(hadTauDecayProducts_.begin(), hadTauDecayProducts_.end(), sortMeasuredHadTauDecayProducts());

    double max_pt = -1.;
    for ( const MeasuredHadTauDecayProduct& hadTauDecayProduct : hadTauDecayProducts_ )
    {
      if ( hadTauDecayProduct.charge() != 0 && hadTauDecayProduct.pt() > max_pt )
      {
        leadChargedHadron_ = &hadTauDecayProduct;
        max_pt = hadTauDecayProduct.pt();
      }
    }
    if ( !leadChargedHadron_ )
    {
      std::cerr << "ERROR: Failed to find leading charged hadron !!" << std::endl;
      assert(0);
    }
  }
}

void
MeasuredTauLepton::initialize()
{
  // CV: pre-compute frequently accessed quantities to reduce computing time;
  //     relations between pT and p, energy taken from http://en.wikipedia.org/wiki/Pseudorapidity
  p_  = pt_*cosh(eta_);
  px_ = pt_*cos(phi_);
  py_ = pt_*sin(phi_);
  pz_ = pt_*sinh(eta_);
  energy_ = std::sqrt(p_*p_ + preciseVisMass_*preciseVisMass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);
  p3_ = Vector(px_, py_, pz_);

  isLeptonicTauDecay_ = (type_ == MeasuredTauLepton::kTauToElecDecay || type_ == MeasuredTauLepton::kTauToMuDecay);
  isHadronicTauDecay_ = (type_ == MeasuredTauLepton::kTauToHadDecay);
  isPrompt_ = (type_ == MeasuredTauLepton::kPrompt);
}

bool
MeasuredTauLepton::isLeptonicTauDecay() const
{
  return isLeptonicTauDecay_;
}

bool
MeasuredTauLepton::isHadronicTauDecay() const
{
  return isHadronicTauDecay_;
}

bool
MeasuredTauLepton::isPrompt() const
{
  return isPrompt_;
}

//---------------------------------------------------------------------------------------------------
// auxiliary function for printing MeasuredTauLepton objects (for debugging purposes)
std::ostream&
operator<<(std::ostream& os, const std::vector<MeasuredTauLepton>& measuredTauLeptons)
{
  for ( size_t idx = 0; idx < measuredTauLeptons.size(); ++idx )
  {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons[idx];
    os << "measuredTauLepton #" << idx 
       << " (type = " << measuredTauLepton.type_string() << ", charge = " << measuredTauLepton.charge() << "):" 
       << " Pt = " << measuredTauLepton.pt() << ","
       << " eta = " << measuredTauLepton.eta() << " (theta = " << measuredTauLepton.p3().theta() << ")" << "," 
       << " phi = " << measuredTauLepton.phi() << ","
       << " mass = " << measuredTauLepton.mass() << std::endl;
  }
  return os;
}
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------
// auxiliary class for sorting MeasuredTauLepton objects
bool
sortMeasuredTauLeptons::operator() (const MeasuredTauLepton& measuredTauLepton1, const MeasuredTauLepton& measuredTauLepton2)
{
  // sort tau decay products into "leg1" (first daughter) and "leg2" (second daughter)
  // for the choice of "leg1", give preference (in order of decreasing priority) to:
  //  - electrons and muons directly originating from lepton-flavor-violating Higgs boson decay
  //  - leptonic over hadronic tau decays 
  //  - tau decay products of higher pT (in case taus decay either both leptonically or both hadronically)
  if ( measuredTauLepton1.isPrompt() && !measuredTauLepton2.isPrompt() ) return true;
  if ( measuredTauLepton2.isPrompt() && !measuredTauLepton1.isPrompt() ) return false;
  // give preference to leptonic tau decays for "leg1";
  // in case taus decay either both leptonically or both hadronically, 
  // give preference to tau decay products of higher pT for "leg1"
  if ( measuredTauLepton1.isLeptonicTauDecay() && measuredTauLepton2.isHadronicTauDecay() ) return true;
  if ( measuredTauLepton2.isLeptonicTauDecay() && measuredTauLepton1.isHadronicTauDecay() ) return false;
  return measuredTauLepton1.pt() > measuredTauLepton2.pt();
}
//---------------------------------------------------------------------------------------------------

}
