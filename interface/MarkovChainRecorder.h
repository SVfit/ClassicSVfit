#ifndef TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h
#define TauAnalysis_ClassicSVfit_svFitHistogramAdapter_h

#include <Math/Functor.h>

namespace classic_svFit
{
  class MarkovChainRecorder : public ROOT::Math::Functor
  {
   public:
    MarkovChainRecorder(unsigned int numDimensions);
    ~MarkovChainRecorder();

    unsigned int
    getNumPoints();

    const double*
    getPoint(unsigned int iPoint);

    double
    getValue(unsigned int iPoint);

   private:
    double
    DoEval(const double* x) const;

    std::vector<std::vector<double>> points_;
    std::vector<double> values_;

    double* x_;
  };
}

#endif
