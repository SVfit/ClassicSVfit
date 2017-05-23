#include "TauAnalysis/ClassicSVfit/interface/SVfitCUBAIntegrator.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <assert.h>

#include "cuba.h"

using namespace classic_svFit;


//
//-------------------------------------------------------------------------------
//
SVfitCUBAIntegrator::SVfitCUBAIntegrator()
  : integrand_(0),
    x_(0){ }
//
//-------------------------------------------------------------------------------
//
SVfitCUBAIntegrator::~SVfitCUBAIntegrator()
{ delete [] x_; }
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::setIntegrand(gPtr_C g, const double* xl, const double* xu, unsigned d){
  numDimensions_ = d;

  delete [] x_;
  x_ = new double[numDimensions_];

  xMin_.resize(numDimensions_);
  xMax_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
  }

  integrand_ = g;
}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::integrate(gPtr_C g, const double* xl, const double* xu, unsigned d, double& integral, double& integralErr)
{
  setIntegrand(g, xl, xu, d);

  if ( !integrand_ ) {
    std::cerr << "<SVfitCUBAIntegrator>:"
              << "No integrand function has been set yet --> ABORTING !!\n";
    assert(0);
  }

  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
    if ( verbosity_ >= 1 ) {
      std::cout << "dimension #" << iDimension << ": min = " << xMin_[iDimension] << ", max = " << xMax_[iDimension] << std::endl;
    }
  }

const int ndim =  numDimensions_;
const int ncomp = 1;
void *userdata = NULL;
const int nvec = 1;
const cubareal epsrel = 1E-3;
const cubareal epsabs = 1E-12;
const int flags = 2;
const int seed = 0;
const int mineval = 0;
const int maxeval = 50000;
const int nstart =  1000;
const int nincrease = 500;
const int nbatch = 1000;
const int gridno =  0;
const char *statefile = NULL;
void *spin = NULL;
int neval;
int fail;
double probV[10];
double integralV[10];
double integralErrV[10];


Vegas(ndim, ncomp,
    integrand_,
    userdata, nvec,
    epsrel, epsabs,
    flags,  seed,
    mineval, maxeval,
    nstart, nincrease, nbatch,
    gridno, statefile, spin,
    &neval, &fail,
    integralV, integralErrV, probV);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for(unsigned int comp = 0; comp < ncomp; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integralV[comp], integralErrV[comp], probV[comp]);

}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::print(std::ostream& stream) const{
  stream<<"VfitCUBAIntegrator::print"<<std::endl;}
//
//-------------------------------------------------------------------------------
//
