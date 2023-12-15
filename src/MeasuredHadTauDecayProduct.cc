#include "TauAnalysis/ClassicSVfit/interface/MeasuredHadTauDecayProduct.h"

#include <TMath.h>

namespace classic_svFit
{

MeasuredHadTauDecayProduct::MeasuredHadTauDecayProduct()
  : charge_(0)
  , pt_(0.)
  , eta_(0.)
  , phi_(0.)
  , mass_(0.)
{
  initialize();
}

MeasuredHadTauDecayProduct::MeasuredHadTauDecayProduct(int charge, double pt, double eta, double phi, double mass)
  : charge_(charge)
  , pt_(pt)
  , eta_(eta)
  , phi_(phi)
  , mass_(mass)
{
  initialize();
}

MeasuredHadTauDecayProduct::MeasuredHadTauDecayProduct(const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct)
  : charge_(measuredHadTauDecayProduct.charge())
  , pt_(measuredHadTauDecayProduct.pt())
  , eta_(measuredHadTauDecayProduct.eta())
  , phi_(measuredHadTauDecayProduct.phi())
  , mass_(measuredHadTauDecayProduct.mass())
{
  initialize();
}

MeasuredHadTauDecayProduct::~MeasuredHadTauDecayProduct()
{}

int MeasuredHadTauDecayProduct::charge() const 
{ 
  return charge_; 
}

double MeasuredHadTauDecayProduct::pt() const 
{ 
  return pt_; 
}

double MeasuredHadTauDecayProduct::eta() const 
{ 
  return eta_; 
}

double MeasuredHadTauDecayProduct::phi() const 
{ 
  return phi_; 
}

double MeasuredHadTauDecayProduct::mass() const 
{ 
  return mass_; 
}

double MeasuredHadTauDecayProduct::energy() const 
{ 
  return energy_; 
}

double MeasuredHadTauDecayProduct::px() const 
{ 
  return px_; 
}

double MeasuredHadTauDecayProduct::py() const 
{ 
  return py_; 
}

double MeasuredHadTauDecayProduct::pz() const 
{ 
  return pz_; 
}

double MeasuredHadTauDecayProduct::p() const 
{
  return p_; 
}

LorentzVector MeasuredHadTauDecayProduct::p4() const 
{
  return p4_; 
}

Vector MeasuredHadTauDecayProduct::p3() const 
{ 
  return p3_; 
}

void MeasuredHadTauDecayProduct::initialize()
{
  // CV: relations between pT and p, energy taken from http://en.wikipedia.org/wiki/Pseudorapidity
  p_  = pt_*TMath::CosH(eta_);
  px_ = pt_*TMath::Cos(phi_);
  py_ = pt_*TMath::Sin(phi_);
  pz_ = pt_*TMath::SinH(eta_);
  energy_ = TMath::Sqrt(p_*p_ + mass_*mass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);
  p3_ = Vector(px_, py_, pz_);
}

//---------------------------------------------------------------------------------------------------
// auxiliary function for printing MeasuredTauLepton objects (for debugging purposes)
std::ostream&
operator<<(std::ostream& os, const std::vector<MeasuredHadTauDecayProduct>& measuredHadTauDecayProducts)
{
  for ( size_t idx = 0; idx < measuredHadTauDecayProducts.size(); ++idx )
  {
    const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct = measuredHadTauDecayProducts[idx];
    os << "measuredHadTauDecayProduct #" << idx << ":" 
       << " (charge = " << measuredHadTauDecayProduct.charge() << "):" 
       << " Pt = " << measuredHadTauDecayProduct.pt() << ","
       << " eta = " << measuredHadTauDecayProduct.eta() << " (theta = " << measuredHadTauDecayProduct.p3().theta() << ")" << "," 
       << " phi = " << measuredHadTauDecayProduct.phi() << ","
       << " mass = " << measuredHadTauDecayProduct.mass() << std::endl;
  }
  return os;
}
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------
// auxiliary class for sorting MeasuredHadTauDecayProduct objects
bool sortMeasuredHadTauDecayProducts::operator() (const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct1, const MeasuredHadTauDecayProduct& measuredHadTauDecayProduct2)
{
  return ( measuredHadTauDecayProduct1.pt() > measuredHadTauDecayProduct2.pt() );
}
//---------------------------------------------------------------------------------------------------

}
