#include "TauAnalysis/ClassicSVfit/interface/comp_PCA_line2point.h"

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

Point
comp_PCA_line2point(const Point& p1, const Vector& v1,
                    const Point& p2,
                    const TMatrixD& covInv,
                    double lambdaMin, double lambdaMax)
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
