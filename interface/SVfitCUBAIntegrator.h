
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

#include "cuba.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

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

void integrate(integrand_t g, const double* xl, const double* xu, unsigned d, double& integral, double& integralErr);

void integrateVector(integrand_t g, const double* xl, const double* xu, unsigned d,
                     double integral[],
                     double integralErr[]);

void setNumberOfComponents(unsigned int);

void print(std::ostream&) const;

protected:

void setIntegrand(integrand_t, const double*, const double*, unsigned);

int verbosity_;     // flag to enable/disable debug output

integrand_t integrand_;

/// parameters defining integration region
///  xMin:          lower boundaries of integration region
///  xMax:          upper boundaries of integration region
double* x_;
std::vector<double> xMin_;     // index = dimension
std::vector<double> xMax_;     // index = dimension


///Cuba library configuration parameters.

int verbose = 0;

///the number of dimensions of the integral.
int ndim;

///the number of components of the integrand.
int ncomp = 1;

///user data passed to the integrand.
void *userdata = NULL;

///the maximum number of points to be given to the integrand routine in each invocation.
int nvec = 1;

///the requested relative and absolute accuracies.
///The integrator tries to find an estimate Iˆ for the integral I which for every
//component c fulfills |Iˆc − Ic |< max(epsabs, epsrel*|Ic |).
cubareal epsrel = 1E-3;
cubareal epsabs = 1E-12;

//flags governing the integration.
int flags = 2;

///the seed for the pseudo-random-number generator.
///seed = 0 corresponds to Sobol (quasi-random)
int seed = 0;

//the minimum number of integrand evaluations required
int mineval = 0;

///the (approximate) maximum number of integrand evaluations allowed
int maxeval = 2000;

//a filename for storing the internal state.
char *statefile = NULL;

///the ‘spinning cores’ pointer. Used for pararell computations
void *spin = NULL;

///VEGAS specific parameters
/// the number of integrand evaluations per iteration to startwith.
int nstart =  500;

///the increase in the number of integrand evaluations per iteration
int nincrease = 500;

///the batch size for sampling.
int nbatch = 1000;

///the slot in the internal grid table
int gridno =  0;
//////////////////////

///Suave specific parameters
///the number of new integrand evaluations in each subdivision
int nnew=1000;

///the minimum number of samples a former pass must contribute to a subregion
///to be considered in that region’s compound integral value.
int nmin=2;

///the parameter p in Eq. (1), od Cuba manual
double flatness=25;
//////////////////////

///Divonne specific parameters
/// determines sampling in the partitioning phase. Value taken from Cuba example.
int key1=47;

///determines sampling in the final integration phase
int key2=1;

///sets the strategy for the refinement phase
int key3=1;

///controls the thoroughness of the partitioning phase
int maxpass=5;

//the width of the border of the integration region
double border=0.;

///the maximum χ2 value a single subregion is allowed to have in the final integration phase.
double maxchisq=10.;

///a bound, given as the fraction of the requested error of the entire integral,
///which determines whether it is worthwhile further examining a region that failed the Chi2 test.
double mindeviation=.25;

///a list of points where the integrand might have peaks
double *xgiven = NULL;

///the number of points in the xgiven array
int ngiven=0;

///the peak-finder subroutine
peakfinder_t peakfinder = NULL;

///the leading dimension of xgiven, i.e. the offset between one point and the next in memory.
int ldxdiven=ndim;

///the maximum number of extra points the peak-finder subroutine will return.
///If nextra is zero, peakfinder is not called
int nextra=0;
//////////////////////

///Cuhre specific parameters
///chooses the basic integration rule
int key = 0;
//////////////////////

///Interagtion output variables.
///Vectors holdig integration result
double integralResult[classic_svFit::maxNumberOfComponents];
double integralError[classic_svFit::maxNumberOfComponents];
double integralProbability[classic_svFit::maxNumberOfComponents];

//the actual number of subregions needed (not present in Vegas).
int nregions = 0;

///the actual number of integrand evaluations needed.
int neval = 0;

///an error flag:
int fail = -1;

};
}

#endif
