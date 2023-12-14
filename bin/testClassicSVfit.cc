
/**
   \class testClassicSVfit testClassicSVfit.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredEvent.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredMEt.h"
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
  double measuredMEtPx =  11.7491;
  double measuredMEtPy = -51.9172;

  // define MET covariance
  TMatrixD covMEt(2, 2);
  covMEt[0][0] =  787.352;
  covMEt[1][0] = -178.63;
  covMEt[0][1] = -178.63;
  covMEt[1][1] =  179.545;

  // define lepton four vectors
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToElecDecay, +1, 33.7393, 0.9409,  -0.541458, 0.51100e-3)); // tau -> electron decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  -1, 25.7322, 0.618228, 2.79362,  0.13957, 0)); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
                    10 three-prong without neutral pions
  */

  MeasuredMEt measuredMEt(measuredMEtPx, measuredMEtPy, covMEt);
  MeasuredEvent measuredEvent(measuredTauLeptons, { measuredMEt });

  int verbosity = 1;
  ClassicSVfit svFitAlgo(verbosity);
#ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo.setHadTauTF(hadTauTF);
  //svFitAlgo.enableHadTauTF();
#endif

  svFitAlgo.enableLogM(6.);
  //svFitAlgo.disableLogM();
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  svFitAlgo.integrate(measuredEvent);
  bool isValidSolution_1stRun = svFitAlgo.isValidSolution();
  double mass_1stRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr_1stRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass_1stRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr_1stRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution_1stRun ) {
    std::cout << "found valid solution: mass = " << mass_1stRun << " +/- " << massErr_1stRun << " (expected value = 115.746 +/- 87.0011),"
              << " transverse mass = " << transverseMass_1stRun << " +/- " << transverseMassErr_1stRun << " (expected value = 114.242 +/- 85.8296)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass_1stRun - 115.746) / 115.746) > 0.001) return 1;
  if (std::abs((massErr_1stRun - 87.001) / 87.0011) > 0.001) return 1;
  if (std::abs((transverseMass_1stRun - 114.242) / 114.242) > 0.001) return 1;
  if (std::abs((transverseMassErr_1stRun - 85.8296) / 85.8296) > 0.001) return 1;
 
  // re-run with mass constraint
  double massContraint = 125.06;
  std::cout << "\n\nTesting integration with ditau mass constraint set to " << massContraint << std::endl;
  svFitAlgo.setLikelihoodFileName("testClassicSVfit_withMassContraint.root");
  svFitAlgo.setDiTauMassConstraint(massContraint);
  svFitAlgo.integrate(measuredEvent);
  bool isValidSolution_2ndRun = svFitAlgo.isValidSolution();
  double mass_2ndRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr_2ndRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass_2ndRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr_2ndRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution_2ndRun ) {
    std::cout << "found valid solution: mass = " << mass_2ndRun << " +/- " << massErr_2ndRun << " (expected value = 124.646 +/- 1.27575),"
              << " transverse mass = " << transverseMass_2ndRun << " +/- " << transverseMassErr_2ndRun << " (expected value = 123.026 +/- 1.19297)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass_2ndRun - 124.646) / 124.646) > 0.001) return 1;
  if (std::abs((massErr_2ndRun - 1.27575) / 1.27575) > 0.001) return 1;
  if (std::abs((transverseMass_2ndRun - 123.026) / 123.026) > 0.001) return 1;
  if (std::abs((transverseMassErr_2ndRun - 1.19297) / 1.19297) > 0.001) return 1;
  
  std::cout << std::endl;
  std::cout << "*****************************************************************************************************************************************" << std::endl;
  std::cout << "* If you use this code, please cite:                                                                                                    *" << std::endl;
  std::cout << "*   L. Bianchini, B. Calpas, J. Conway, A. Fowlie, L. Marzola, L. Perrini, C. Veelken,                                                  *" << std::endl;
  std::cout << "*   \"Reconstruction of the Higgs mass in events with Higgs bosons decaying into a pair of tau leptons using matrix element techniques\", *" << std::endl;
  std::cout << "*   Nucl. Instrum. Meth. A 862 (2017) 54                                                                                                *" << std::endl;
  std::cout << "*****************************************************************************************************************************************" << std::endl;
  std::cout << std::endl;

  return 0;
}
