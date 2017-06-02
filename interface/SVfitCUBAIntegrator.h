
#ifndef TauAnalysis_ClassicSVfit_SVfitCUBAIntegrator_h
#define TauAnalysis_ClassicSVfit_SVfitCUBAIntegrator_h

/** \class SVfitCUBAIntegrator
 *
 * Wrapper class for CUBA integrators.
 * http://www.feynarts.de/cuba/
 *
 * \author Artur Kalinowski, Uniwersity of Warsaw, Warsaw
 *
 */

#include <Math/Functor.h>

#include <vector>
#include <string>
#include <iostream>

namespace classic_svFit
{
  class SVfitCUBAIntegrator
  {
   public:
    SVfitCUBAIntegrator(unsigned int verbosity);
    ~SVfitCUBAIntegrator();

    /// compute integral of function g
    /// the points xl and xh represent the lower left and upper right corner of a Hypercube in d-dimensional integration space
    typedef int (*integrand_t)(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);

    void integrate(integrand_t g, const double* xl, const double* xu,
      unsigned d, double& integral, double& integralErr, float testMass = -1);

    void print(std::ostream&) const;

  protected:

    void setIntegrand(integrand_t, const double*, const double*, unsigned);

    int verbosity_; // flag to enable/disable debug output

    integrand_t integrand_;

    /// parameters defining integration region
    ///  numDimensions: dimensionality of integration region (Hypercube)
    ///  xMin:          lower boundaries of integration region
    ///  xMax:          upper boundaries of integration region
    ///  initMode:      flag indicating how initial position of Markov Chain is chosen (uniform/Gaus distribution)
    unsigned numDimensions_;
    double* x_;
    std::vector<double> xMin_; // index = dimension
    std::vector<double> xMax_; // index = dimension

  };
}

#endif
