#include "TauAnalysis/ClassicSVfit/interface/comp_PCA.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h" // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include <TVectorD.h>                            // TVectorD

#include <assert.h>                              // assert()

using namespace classic_svFit;

namespace
{
  std::pair<Point, Point>
  comp_PCA_line2line(
    const Point& p1, const Vector& v1,
    const Point& p2, const Vector& v2,
    const TMatrixD& covInv,
    double lambda1Min = -1.e+6, double lambda1Max = +1.e+6, double lambda2Min = -1.e+6, double lambda2Max = +1.e+6)
  {
std::cout << "<comp_PCA_line2line>:" << std::endl;
std::cout << " p1: x = " << p1.x() << ", y = " << p1.y() << ", z = " << p1.z() << std::endl;
std::cout << " v1: theta = " << v1.theta() << ", phi = " << p1.phi() << std::endl;
std::cout << " p2: x = " << p2.x() << ", y = " << p2.y() << ", z = " << p2.z() << std::endl;
std::cout << " v2: theta = " << v2.theta() << ", phi = " << p2.phi() << std::endl;
    // CV: compute point of closest approach (PCA) between two straight lines (p1 + lambda1*v1) and (p2 + lambda2*v2) in three dimensions;
    //     code based on https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines

    TVectorD e1 = convert_to_mathVector(v1.unit());
    TVectorD covInv_times_e1 = covInv*e1;
    TVectorD e2 = convert_to_mathVector(v2.unit());
    TVectorD covInv_times_e2 = covInv*e2;
    TMatrixD A(2,2);
    A(0,0) =   e1*covInv_times_e1;
    A(0,1) = -(e1*covInv_times_e2);
    A(1,0) =   e2*covInv_times_e1;
    A(1,1) = -(e2*covInv_times_e2);
    if ( A.Determinant() == 0. )
    {
      std::cerr << "ERROR: Failed to invert matrix A !!" << std::endl;
      assert(0);
    }
    TMatrixD Ainv(TMatrixD::kInverted, A);

    TVectorD p2_minus_p1 = convert_to_mathVector(p2 - p1);
    TVectorD covInv_times_p2_minus_p1 = covInv*p2_minus_p1;
    TVectorD c(2);
    c(0) = e1*covInv_times_p2_minus_p1;
    c(1) = e2*covInv_times_p2_minus_p1;

    TVectorD lambda = Ainv*c;
    double lambda1 = lambda(0);
    if ( lambda1 < lambda1Min ) lambda1 = lambda1Min;
    if ( lambda1 > lambda1Max ) lambda1 = lambda1Max;
    double lambda2 = lambda(1);
    if ( lambda2 < lambda2Min ) lambda2 = lambda2Min;
    if ( lambda2 > lambda2Max ) lambda2 = lambda2Max;

    Point pca1 = p1 + lambda1*convert_to_recoVector(e1);
    Point pca2 = p2 + lambda2*convert_to_recoVector(e2);
std::cout << "pca1: x = " << pca1.x() << ", y = " << pca1.y() << ", z = " << pca1.z() << std::endl;
std::cout << "pca2: x = " << pca2.x() << ", y = " << pca2.y() << ", z = " << pca2.z() << std::endl;
    return std::pair<Point, Point>(pca1, pca2);
  }

  Point
  comp_PCA_line2point(
    const Point& p1, const Vector& v1,
    const Point& p2,
    const TMatrixD& covInv,
    double lambdaMin = -1.e+6, double lambdaMax = +1.e+6)
  {
    // CV: compute point of closest approach (PCA) between straight line (p1 + lambda*v1) and point (p2) in three dimensions.

    TVectorD e1 = convert_to_mathVector(v1.unit());
    TVectorD covInv_times_e1 = covInv*e1;
  
    TVectorD p2_minus_p1 = convert_to_mathVector(p2 - p1);
    double lambda = (p2_minus_p1*covInv_times_e1)/(e1*covInv_times_e1);
    if ( lambda < lambdaMin ) lambda = lambdaMin;
    if ( lambda > lambdaMax ) lambda = lambdaMax;

    Point pca = p1 + lambda*convert_to_recoVector(e1);
    return pca;
  }
}

Point
classic_svFit::comp_PCA(
  const LorentzVector& tauP4,
  const MeasuredTauLepton& measuredTauLepton, const MeasuredHadTauDecayProduct& measuredLeadChargedHadron,
  const Point& primaryVertex, const Point& decayVertex, const TMatrixD& decayVertexCovInv)
{
  Point pca;
  if ( measuredTauLepton.decayMode() == reco::PFTau::kThreeProng0PiZero || 
       measuredTauLepton.decayMode() == reco::PFTau::kThreeProng1PiZero )
  {
    pca = ::comp_PCA_line2point(
            primaryVertex, tauP4.Vect(),
            decayVertex, decayVertexCovInv,
            0., 1.e+6);
  }
  else
  {
    pca = ::comp_PCA_line2line(
            primaryVertex, tauP4.Vect(),
            decayVertex, measuredLeadChargedHadron.p3(), decayVertexCovInv,
            0., 1.e+6, -1.e+6, +1.e+6).first;
  }
  return pca;
}
