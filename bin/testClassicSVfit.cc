
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

  double mass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution ) {
    std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " (expected value = 115.746 +/- 92.5252),"
              << " transverse mass = " << transverseMass << " +/- " << transverseMassErr << " (expected value = 114.242 +/- 91.2066)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass - 115.746) / 115.746) > 0.001) return 1;
  if (std::abs((massErr - 92.5252) / 92.5252) > 0.001) return 1;
  if (std::abs((transverseMass - 114.242) / 114.242) > 0.001) return 1;
  if (std::abs((transverseMassErr - 91.2066) / 91.2066) > 0.001) return 1;
  
  // re-run with mass constraint
  double massContraint = 125.06;
  std::cout << "\n\nTesting integration with di tau mass constraint set to " << massContraint << std::endl;
  svFitAlgo.setLikelihoodFileName("testClassicSVfit_withMassContraint.root");
  svFitAlgo.setDiTauMassConstraint(massContraint);
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  isValidSolution = svFitAlgo.isValidSolution();
  mass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
  massErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  transverseMass = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  transverseMassErr = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution ) {
    std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " (expected value = 124.646 +/- 1.23027),"
              << " transverse mass = " << transverseMass << " +/- " << transverseMassErr << " (expected value = 123.026 +/- 1.1574)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass - 124.646) / 124.646) > 0.001) return 1;
  if (std::abs((massErr - 1.23027) / 1.23027) > 0.001) return 1;
  if (std::abs((transverseMass - 123.026) / 123.026) > 0.001) return 1;
  if (std::abs((transverseMassErr - 1.1574) / 1.1574) > 0.001) return 1;
  
  // re-run with classic_svFit::TauTauHistogramAdapter
  std::cout << "\n\nTesting integration with classic_svFit::TauTauHistogramAdapter" << std::endl;
  ClassicSVfit svFitAlgo2(verbosity);
  svFitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
  svFitAlgo2.addLogM_fixed(true, 6.);
  svFitAlgo2.setLikelihoodFileName("testClassicSVfit_TauTauHistogramAdapter.root");
  svFitAlgo2.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  isValidSolution = svFitAlgo2.isValidSolution();
  classic_svFit::LorentzVector tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svFitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
  classic_svFit::LorentzVector tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svFitAlgo2.getHistogramAdapter())->GetFittedTau2LV();

  if ( isValidSolution ) {
    std::cout << "found valid solution: pT(tau1) = " << tau1P4.Pt() << " (expected value = 102.508),"
              << "                      pT(tau2) = " << tau2P4.Pt() << " (expected value = 27.019)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((tau1P4.Pt() - 102.508) / 102.508) > 0.001) return 1;
  if (std::abs((tau2P4.Pt() - 27.019) / 27.019) > 0.001) return 1;

  return 0;
}
