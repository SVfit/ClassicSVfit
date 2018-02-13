#include "TauAnalysis/ClassicSVfit/interface/SVfitCUBAIntegrator.h"

#include <TMath.h>
#include <TUUID.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <assert.h>

using namespace classic_svFit;

//
//-------------------------------------------------------------------------------
//
SVfitCUBAIntegrator::SVfitCUBAIntegrator(unsigned int verbosity, unsigned int maxObjFunctionCalls)
        : integrand_(0), verbosity_(verbosity), x_(0){

        ///Disable multithread calculation, as this can affect the grid running.
        cubacores(0,1);

        ///Set approximate maximum number of integrand calls
        maxeval = maxObjFunctionCalls;

}
//
//-------------------------------------------------------------------------------
//
SVfitCUBAIntegrator::~SVfitCUBAIntegrator()
{
        delete [] x_;
}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::setIntegrand(integrand_t g, const double* xl, const double* xu, unsigned d){

        ndim = d;

        delete [] x_;
        x_ = new double[ndim];

        xMin_.resize(ndim);
        xMax_.resize(ndim);
        for ( unsigned iDimension = 0; iDimension < ndim; ++iDimension ) {
                xMin_[iDimension] = xl[iDimension];
                xMax_[iDimension] = xu[iDimension];
                if ( verbosity_ >= 2 ) {
                        std::cout << "dimension #" << iDimension << ": min = " << xMin_[iDimension] << ", max = " << xMax_[iDimension] << std::endl;
                }
        }

        integrand_ = g;

        if ( !integrand_ ) {
                std::cerr << "<SVfitCUBAIntegrator>:"
                          << "No integrand function has been set yet --> ABORTING !!\n";
                assert(0);
        }
}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::setNumberOfComponents(unsigned int aNumberOfComp){

        ncomp = aNumberOfComp;

}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::integrate(integrand_t g, const double* xl, const double* xu,
                                    unsigned d, double& integral, double& integralErr){

        setIntegrand(g, xl, xu, d);

        Divonne(ndim, ncomp, integrand_, userdata, nvec,
                epsrel, epsabs, verbose, seed,
                mineval, maxeval, key1, key2, key3, maxpass,
                border, maxchisq, mindeviation,
                ngiven, ldxdiven, xgiven, nextra, peakfinder,
                statefile, spin,
                &nregions, &neval, &fail,
                integralResult, integralError, integralProbability);

/*
   Vegas(ndim, ncomp, integrand_,
    userdata, nvec, epsrel, epsabs,
    verbose,  seed, mineval, maxeval,
    nstart, nincrease, nbatch,
    gridno, statefile, spin,
    &neval, &fail,
    integralV, integralErrV, probV);

   Suave(ndim, ncomp, integrand_, userdata, nvec,
    epsrel, epsabs, verbose, seed,
    mineval, maxeval, nnew, nmin, flatness,
    statefile, spin,
    &nregions, &neval, &fail, integralResult, integralError, integralProbability);

    Cuhre(ndim, ncomp, integrand_, userdata, nvec,
    epsrel, epsabs, verbose,
    mineval, maxeval, key,
    statefile, spin,
    &nregions, &neval, &fail, integralResult, integralError, integralProbability);
 */

        integral = integralResult[0];
        integralErr = integralError[0];
}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::integrateVector(integrand_t g, const double* xl, const double* xu, unsigned d,
                                          double aIntegral[],
                                          double aIntegralErr[]){

        setIntegrand(g, xl, xu, d);

        Divonne(ndim, ncomp, integrand_, userdata, nvec,
                epsrel, epsabs, verbose, seed,
                mineval, maxeval, key1, key2, key3, maxpass,
                border, maxchisq, mindeviation,
                ngiven, ldxdiven, xgiven, nextra, peakfinder,
                statefile, spin,
                &nregions, &neval, &fail,
                aIntegral, aIntegralErr, integralProbability);

}
//
//-------------------------------------------------------------------------------
//
void SVfitCUBAIntegrator::print(std::ostream& stream) const {
        stream<<"SVfitCUBAIntegrator::print"<<std::endl;
}
//
//-------------------------------------------------------------------------------
//
