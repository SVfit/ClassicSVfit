
/**
   \class testClassicSVfit testClassicSVfit.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

using namespace classic_svFit;

int main(int argc, char* argv[])
{
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
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution = svFitAlgo.isValidSolution();
  double mass = svFitAlgo.mass();
  double massErr = svFitAlgo.massErr();
  double transverseMass = svFitAlgo.transverseMass();
  double transverseMassErr = svFitAlgo.transverseMassErr();
  if ( isValidSolution ) {
    std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " (expected value = 115.746 +/- 88.6066),"
              << " transverse mass = " << transverseMass << " +/- " << transverseMassErr << " (expected value = 114.242 +/- 87.4277)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }

  return 0;
}

