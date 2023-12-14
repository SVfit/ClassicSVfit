#ifndef TauAnalysis_ClassicSVfit_comp_PCA_line2line_h
#define TauAnalysis_ClassicSVfit_comp_PCA_line2line_h

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // Point, Vector

#include <TMatrixD.h>                                             // TMatrixD

#include <utility>                                                // std::pair<>

namespace classic_svFit
{
  std::pair<Point, Point>
  comp_PCA_line2line(const Point& p1, const Vector& v1,
                     const Point& p2, const Vector& v2,
                     const TMatrixD& covInv,
                     double lambda1Min = -1.e+6, double lambda1Max = +1.e+6, double lambda2Min = -1.e+6, double lambda2Max = +1.e+6);
}

#endif


