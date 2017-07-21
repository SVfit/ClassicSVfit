
/**
   \class testClassicSVfit testClassicSVfit.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "TH1F.h"

using namespace classic_svFit;

typedef double v8si __attribute__ ((vector_size (8*sizeof(double))));

v8si xMin = {0,0,0,0,0,0,0,0};
v8si xMax = {1,1,1,1,1,1,1,1};
v8si x = {0,0,0,0,0,0,0,0};

std::vector<double> xMin_ = {0,0,0,0,0,0,0,0};
std::vector<double> xMax_ = {1,1,1,1,1,1,1,1};
std::vector<double> x_ = {0,0,0,0,0,0,0,0};

void updateXSIMD(const v8si & q){

x = xMin*(1-q) + xMax;

}

void updateX(const std::vector<double>& q)
{
  for ( unsigned iDimension = 0; iDimension < 6; ++iDimension ) {
    const double & q_i = q[iDimension];
    x_[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }
}

int main(int argc, char* argv[])
{
/*
v8si q = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
std::vector<double> q_ = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

for(unsigned int i=0;i<1E6;++i){
updateXSIMD(q);
updateX(q_);
}
return 0;
*/

  /*
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMETx =  11.7491;
  double measuredMETy = -51.9172;

  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] =  787.352;
  covMET[1][0] = -178.63;
  covMET[0][1] = -178.63;
  covMET[1][1] =  179.545;

  // define lepton four vectors
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, 33.7393, 0.9409,  -0.541458, 0.51100e-3)); // tau -> electron decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  25.7322, 0.618228, 2.79362,  0.13957, 0)); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
        10 three-prong without neutral pions
  */

  int verbosity = 1;
  ClassicSVfit svFitAlgo(verbosity);
#ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo.setHadTauTF(hadTauTF);
  //svFitAlgo.enableHadTauTF();
#endif
  //svFitAlgo.addLogM_fixed(false);
  svFitAlgo.addLogM_fixed(true, 6.);
  //svFitAlgo.addLogM_dynamic(true, "(m/1000.)*15.");
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");

  TFile testFile("Test.root","RECREATE");
  TH1F *hCubaMass = new TH1F("hCubaMass","",30,100,130);
  TH1F *hMCMass = new TH1F("hMCMass","",30,100,130);

  for(unsigned int iTry=0;iTry<1;++iTry){
    svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
    double mass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
    hMCMass->Fill(mass);
  }
  bool isValidSolution = svFitAlgo.isValidSolution();

  double lMax = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMassLmax();
  double mass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  float maxMass = 0;
  for(unsigned int iTry=0;iTry<1;++iTry){
     maxMass = svFitAlgo.integrateCuba(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
     hCubaMass->Fill(maxMass);
   }
   testFile.Write();

  if ( isValidSolution ) {
    std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " (expected value = 115.746 +/- 88.6115),"
              << " transverse mass = " << transverseMass << " +/- " << transverseMassErr << " (expected value = 114.242 +/- 87.4328)"
              << std::endl<<"  solution with Cuba: mass = "<<maxMass << " (expected value = 110.169)"
              << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }

  return 0;
}
