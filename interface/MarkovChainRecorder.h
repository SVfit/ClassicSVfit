#ifndef TauAnalysis_ClassicSVfit_MarkovChainRecorder_h
#define TauAnalysis_ClassicSVfit_MarkovChainRecorder_h

#include <Math/Functor.h> // ROOT::Math::Functor

#include <string>         // std::string
#include <vector>         // std::vector<>

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

    unsigned int numDimensions_;

    mutable std::vector<std::vector<double>> points_;
    mutable std::vector<double> values_;

    std::vector<double> x_;
  };
}

#endif
