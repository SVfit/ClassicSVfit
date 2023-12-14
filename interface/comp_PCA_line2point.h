#ifndef TauAnalysis_ClassicSVfit_comp_PCA_line2point_h
#define TauAnalysis_ClassicSVfit_comp_PCA_line2point_h

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // Point, Vector

#include <TMatrixD.h>                                             // TMatrixD

namespace classic_svFit
{
  Point
  comp_PCA_line2point(const Point& p1, const Vector& v1,
                      const Point& p2,
                      const TMatrixD& covInv,
                      double lambdaMin = -1.e+6, double lambdaMax = +1.e+6);
}

#endif

