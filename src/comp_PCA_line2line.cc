#include "TauAnalysis/ClassicSVfit/interface/comp_PCA_line2line.h"

#include <TVectorD.h> // TVectorD

#include <assert.h>   // assert()

using namespace classic_svFit;

namespace
{
  template <typename T>
  TVectorD
  convert_to_mathVector(const T& v)
  {
    TVectorD retVal(3);
    retVal(0) = v.x();
    retVal(1) = v.y();
    retVal(2) = v.z();
    return retVal;
  }

  Vector
  convert_to_recoVector(const TVectorD& v)
  {
    return Vector(v(0), v(1), v(2));
  }
}

namespace classic_svFit
{

std::pair<Point, Point>
comp_PCA_line2line(const Point& p1, const Vector& v1,
                   const Point& p2, const Vector& v2,
                   const TMatrixD& covInv,
                   double lambda1Min, double lambda1Max, double lambda2Min, double lambda2Max)
{
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
  return std::pair<Point, Point>(pca1, pca2);
}

}
