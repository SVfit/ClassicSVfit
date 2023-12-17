#ifndef TauAnalysis_ClassicSVfit_comp_PCA_h
#define TauAnalysis_ClassicSVfit_comp_PCA_h

#include "TauAnalysis/ClassicSVfit/interface/MeasuredHadTauDecayProduct.h" // MeasuredHadTauDecayProduct
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"          // MeasuredTauLepton
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"          // Point, Vector

#include <TMatrixD.h>                                                      // TMatrixD

namespace classic_svFit
{
  Point
  comp_PCA(const LorentzVector& tauP4,
           const MeasuredTauLepton& measuredTauLepton, const MeasuredHadTauDecayProduct& measuredLeadChargedHadron,
           const Point& primaryVertex, const Point& decayVertex, const TMatrixD& decayVertexCovInv);
}

#endif
