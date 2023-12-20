#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

#include <TMath.h>    // TMath::Floor(), TMath::Nint()

#include <cmath>      // log10(), pow(), std::fabs(), std::sqrt()

using namespace classic_svFit;

double
classic_svFit::roundToNdigits(double x, int n)
{
  double tmp = pow(10., n);
  if ( x != 0. )
  {
    tmp /= pow(10., TMath::Floor(log10(std::fabs(x))));
  }
  double x_rounded = TMath::Nint(x*tmp)/tmp;
  return x_rounded;
}

TMatrixD
classic_svFit::roundToNdigits(const TMatrixD& m, int n)
{
  int nRows = m.GetNrows();
  int nColumns = m.GetNcols();
  TMatrixD m_rounded(nRows, nColumns);
  for ( int iRow = 0; iRow < nRows; ++iRow )
  {
    for ( int iColumn = 0; iColumn < nColumns; ++iColumn )
    {
      m_rounded(iRow,iColumn) = roundToNdigits(m(iRow,iColumn), n);
    }
  }
  return m_rounded;
}

Vector
classic_svFit::normalize(const Vector& p)
{
  double p_x = p.x();
  double p_y = p.y();
  double p_z = p.z();
  double mag2 = square(p_x) + square(p_y) + square(p_z);
  if ( mag2 <= 0. ) return p;
  double mag = std::sqrt(mag2);
  return Vector(p_x/mag, p_y/mag, p_z/mag);
}

double
classic_svFit::compScalarProduct(const Vector& p1, const Vector& p2)
{
  return (p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
}

Vector
classic_svFit::compCrossProduct(const Vector& p1, const Vector& p2)
{
  double p3_x = p1.y()*p2.z() - p1.z()*p2.y();
  double p3_y = p1.z()*p2.x() - p1.x()*p2.z();
  double p3_z = p1.x()*p2.y() - p1.y()*p2.x();
  return Vector(p3_x, p3_y, p3_z);
}

double
classic_svFit::compCosThetaNuNu(double visEn, double visP, double visMass2, double nunuEn, double nunuP, double nunuMass2)
{
  double cosThetaNuNu = (visEn*nunuEn - 0.5*(tauLeptonMass2 - (visMass2 + nunuMass2)))/(visP*nunuP);
  return cosThetaNuNu;
}

double
classic_svFit::compPSfactor_tauToLepDecay(double x, double visEn, double visP, double visMass, double nunuEn, double nunuP, double nunuMass)
{
  double visMass2 = square(visMass);
  double nunuMass2 = square(nunuMass);
  if ( x >= (visMass2/tauLeptonMass2) && x <= 1. && nunuMass2 < ((1. - x)*tauLeptonMass2) ) 
  { 
    // physical solution
    double tauEn_rf = (tauLeptonMass2 + nunuMass2 - visMass2)/(2.*nunuMass);
    double visEn_rf = tauEn_rf - nunuMass;
    if ( !(tauEn_rf >= tauLeptonMass && visEn_rf >= visMass) ) return 0.;
    double I = nunuMass2*(2.*tauEn_rf*visEn_rf - (2./3.)*std::sqrt((square(tauEn_rf) - tauLeptonMass2)*(square(visEn_rf) - visMass2)));
    #ifdef XSECTION_NORMALIZATION
    I *= GFfactor;    
    #endif
    double cosThetaNuNu = classic_svFit::compCosThetaNuNu(visEn, visP, visMass2, nunuEn, nunuP, nunuMass2);
    if ( !(cosThetaNuNu >= (-1. + epsilon) && cosThetaNuNu <= +1.) ) return 0.;
    double PSfactor = (visEn + nunuEn)*I/(8.*visP*square(x)*std::sqrt(square(visP) + square(nunuP) + 2.*visP*nunuP*cosThetaNuNu + tauLeptonMass2));
    //-------------------------------------------------------------------------
    // CV: fudge factor to reproduce literature value for cross-section times branching fraction
    #ifdef XSECTION_NORMALIZATION
    PSfactor *= 2.;
    #endif
    //-------------------------------------------------------------------------
    return PSfactor;
  } 
  else 
  {
    return 0.;
  }
}

double
classic_svFit::compPSfactor_tauToHadDecay(double x, double visEn, double visP, double visMass, double nuEn, double nuP)
{
  double visMass2 = square(visMass);
  if ( x >= (visMass2/tauLeptonMass2) && x <= 1. ) { // physical solution
    double cosThetaNu = classic_svFit::compCosThetaNuNu(visEn, visP, visMass2, nuEn, nuP, 0.);
    if ( !(cosThetaNu >= (-1. + epsilon) && cosThetaNu <= +1.) ) return 0.;
    double PSfactor = (visEn + nuEn)/(8.*visP*square(x)*std::sqrt(square(visP) + square(nuP) + 2.*visP*nuP*cosThetaNu + tauLeptonMass2));
    PSfactor *= 1.0/(tauLeptonMass2 - visMass2);
    //-------------------------------------------------------------------------
    // CV: multiply by constant matrix element,
    //     chosen such that the branching fraction of the tau to decay into hadrons is reproduced
    //const double M2 = 16.*TMath::Pi()*cube(tauLeptonMass)*GammaTauToHad/(tauLeptonMass2 - visMass2);
    //Remove multiplication as it add to execution time, and does not alter the result.
    #ifdef XSECTION_NORMALIZATION
    PSfactor *= M2;
    #endif
    //-------------------------------------------------------------------------
    return PSfactor;
  } 
  else 
  {
    return 0.;
  }
}

namespace classic_svFit
{

integrationParameters::integrationParameters()
{
  reset();
}

integrationParameters::~integrationParameters()
{}

void
integrationParameters::reset()
{
  idx_X_ = -1;
  idx_phi_ = -1;
  idx_VisPtShift_ = -1;
  idx_mNuNu_ = -1;
  idx_flightLength_ = -1;
}

}

namespace
{
  classic_svFit::LorentzVector
  fixMass(const classic_svFit::LorentzVector& p4, double mass)
  {
    double px     = p4.px();
    double py     = p4.py();
    double pz     = p4.pz();
    double energy = std::sqrt(px*px + py*py + pz*pz + mass*mass);
    classic_svFit::LorentzVector p4_fixed(px, py, pz, energy);
    return p4_fixed;
  }
}

LorentzVector
classic_svFit::fixTauMass(const LorentzVector& tauP4)
{
  return fixMass(tauP4, tauLeptonMass);
}

LorentzVector
classic_svFit::fixNuMass(const LorentzVector& nuP4)
{
  return fixMass(nuP4, 0.);
}

std::pair<double,double>
classic_svFit::comp_dmin_and_dmax(const LorentzVector& tauP4, const Vector& flightLength, const TMatrixD& covDecayVertex)
{
  Vector tauP3 = tauP4.Vect();
  TVectorD eTau = convert_to_mathVector(normalize(tauP3));
  double d = std::sqrt(flightLength.mag2());
  double sigma2 = eTau*(covDecayVertex*eTau);
  double gamma = tauP4.energy()/classic_svFit::tauLeptonMass;
  if ( gamma < 1. )
  {
    std::cerr << "WARNING: gamma = " << gamma << " is unphysical, setting it to 1 !!" << std::endl;
    gamma = 1.;
  }
  double gamma_times_cTauLifetime = gamma*classic_svFit::cTauLifetime;
  double dmin = 0.;
  double dmax = 0.;
  if ( sigma2 > 0. )
  {
    double sigma = std::sqrt(sigma2);
    if ( sigma < 1.e-3 ) sigma = 1.e-3;
    dmin = d - 5.*sigma;
    dmax = d + 5.*sigma;
    if ( dmin < 0. )
    {
      dmin = 0.;
      dmax = std::min(d + 5.*sigma, 10.*gamma_times_cTauLifetime);
    }
  }
  else
  {
    dmin = 0.;
    dmax = 10.*gamma_times_cTauLifetime;
  }
  if ( !(dmax > dmin) )
  {
    std::cerr << "ERROR: Failed to compute lower and upper limits on tauFlightLength !!" << std::endl;
    std::cout << "covDecayVertex:" << std::endl;
    covDecayVertex.Print();
    std::cout << "eTau:" << std::endl;
    eTau.Print();
    std::cout << "d = " << d << std::endl;
    std::cout << "gamma*ctau = " << gamma_times_cTauLifetime << std::endl;
    std::cout << "sigma2 = " << sigma2 << std::endl;
    std::cout << "dmin = " << dmin << ", dmax = " << dmax << std::endl;
    assert(0);
  }
  return std::make_pair(dmin, dmax);
}

Vector
classic_svFit::convert_to_recoVector(const TVectorD& v)
{
  return Vector(v(0), v(1), v(2));
}

LorentzVector
classic_svFit::get_beamP4(double beamE, double mBeamParticle)
{
  double beamPx = 0.;
  double beamPy = 0.;
  double beamPz = std::sqrt(square(beamE) - square(mBeamParticle));
  LorentzVector beamP4(beamPx, beamPy, beamPz, beamE);
  return beamP4;
}

Vector
classic_svFit::get_r(const Vector& k, const Vector& h)
{
  double cosTheta = k.Dot(h);
  // CV: allow for small rounding errors
  if ( cosTheta < -1.01 || cosTheta > +1.01 )
  {
    std::cerr << "ERROR: cosTheta = " << cosTheta << " outside physical range !!\n";
    assert(0);
  }
  if ( cosTheta < -1. ) cosTheta = -1.;
  if ( cosTheta > +1. ) cosTheta = +1.;
  double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  Vector r = (h - k*cosTheta)*(1./sinTheta);
  return r;
}

Vector
classic_svFit::get_n(const Vector& k, const Vector& r)
{
  // CV: The ordering of r and k in the cross product has been agreed with Luca on 06/09/2023.
  //     The definition n = r x k has been chosen for consistency with Eq. (2.5) in the paper arXiv:1508.05271,
  //     which Luca and Marco have used in their previous papers on Entanglement.
  //    (Whether one computes the vector n using n = r x k or using n = p x k makes no difference:
  //     in both cases, the vector n refers to the direction perpendicular to the scattering plane
  //     and the vectors { n, r, k } define a right-handed coordinate system)
  Vector n = r.Cross(k);
  return n;
}

void
classic_svFit::get_localCoordinateSystem(const Vector& flightDirection, Vector& r, Vector& n, Vector& k)
{
  k = flightDirection.unit();
  Vector h = get_beamP4().Vect().unit();
  r = get_r(k, h);
  n = get_n(k, r);
}

namespace
{
  TMatrixD
  get_rotationMatrixImp(const Vector& r, const Vector& n, const Vector& k, bool kInverted)
  {
    TMatrixD rotMatrix(3,3);
    std::vector<Vector> rnk = { r, n, k };
    classic_svFit::Vector x(1.,0.,0.);
    classic_svFit::Vector y(0.,1.,0.);
    classic_svFit::Vector z(0.,0.,1.);
    std::vector<Vector> xyz = { x, y, z };
    for ( unsigned int i = 0; i < 3; ++i )
    {
      for ( unsigned int j = 0; j < 3; ++j )
      {
        if ( kInverted )
        {
          const Vector& e_i = xyz[i];
          const Vector& e_j = rnk[j];
          rotMatrix(i,j) = e_i.Dot(e_j);
        }
        else
        {
          const Vector& e_i = rnk[i];
          const Vector& e_j = xyz[j];
          rotMatrix(i,j) = e_i.Dot(e_j);
        }
      }
    }
    return rotMatrix;
  }
}

TMatrixD
classic_svFit::get_rotationMatrix(const Vector& r, const Vector& n, const Vector& k)
{
  // compute rotation matrix for transformation from laboratory frame { x, y, z } to helicity frame { r, n, k }
  return get_rotationMatrixImp(r, n, k, false);
}

TMatrixD
classic_svFit::get_rotationMatrixInv(const Vector& r, const Vector& n, const Vector& k)
{
  // compute rotation matrix for transformation from helicity frame { r, n, k } back to laboratory frame { x, y, z }
  return get_rotationMatrixImp(r, n, k, true);
}

TMatrixD
classic_svFit::rotateCovMatrix(const TMatrixD& cov, const TMatrixD& rotMatrix)
{
  // compute elements of covariance matrix cov in rotated coordinate system, given by the matrix rotMatrix,
  // as described in Example 4.25 in Section 4.9.2 of
  //   V. Blobel and E. Lohrmann "Statistische und numerische Methoden der Datenanalyse".
  TMatrixD rotMatrixT = TMatrixD(TMatrixD::kTransposed, rotMatrix);
  return rotMatrix*cov*rotMatrixT;
}

TMatrixD
classic_svFit::invertMatrix(const std::string& label, const TMatrixD& cov, bool& errorFlag)
{
  TMatrixD covInv;
  errorFlag = false;
  if ( cov.Determinant() == 0. )
  {
    std::cout << label << ":" << std::endl;
    cov.Print();
    std::cerr << "ERROR: Failed to invert matrix cov (det=0) !!" << std::endl;
    errorFlag = true;
    return covInv;
  }
  covInv.ResizeTo(3,3);
  covInv = TMatrixD(TMatrixD::kInverted, cov);
  return covInv;
}

