
/**
   \class testClassicSVfitLFV testClassicSVfitLFV .cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfitLFV .cc"
   \brief Example for computing Higgs boson (H) mass in lepton-flavor-violating H decays, e.g H -> mu tau -> mu tau_h nu (cf. arXiv:1502.07400)
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h" // LorentzVector
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "TH1F.h"
#include "TMath.h"

using namespace classic_svFit;

int main(int argc, char* argv[])
{
  /*
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMEtPx = 17.6851;
  double measuredMEtPy = 23.5161;

  // define MET covariance
  TMatrixD covMEt(2, 2);
  covMEt[0][0] = 284.0;
  covMEt[1][0] =  13.4;
  covMEt[0][1] =  13.4;
  covMEt[1][1] = 255.6;

  //-------------------------------------------------------------------------------------------------
  // CV: this code is a bad example! 
  //     Please store the pT, eta, phi, mass components in your Ntuple directly. 
  //     The conversion from cartesian coordinates to pT, eta, phi, mass may produce non-negligible rounding errors,
  //     in particular on the mass (which is of special importance to SVfit, as it limits the "physical" region in integration space)
  double leg1Px   = -47.5987;
  double leg1Py   = -13.6761;
  double leg1Pz   = -61.7664;
  double leg1Mass =   0.105658;
  double leg1En   = TMath::Sqrt(leg1Px*leg1Px + leg1Py*leg1Py + leg1Pz*leg1Pz + leg1Mass*leg1Mass);
  classic_svFit::LorentzVector leg1P4(leg1Px, leg1Py, leg1Pz, leg1En); 
  
  double leg2Px   =  35.0206;
  double leg2Py   =   9.57334;
  double leg2Pz   =   9.49413;
  double leg2Mass =   1.00231;
  double leg2En   = TMath::Sqrt(leg2Px*leg2Px + leg2Py*leg2Py + leg2Pz*leg2Pz + leg2Mass*leg2Mass);
  classic_svFit::LorentzVector leg2P4(leg2Px, leg2Py, leg2Pz, leg2En); 
  //-------------------------------------------------------------------------------------------------

  // define lepton four vectors
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kPrompt, +1, leg1P4.pt(), leg1P4.eta(), leg1P4.phi(), leg1P4.mass())); // "prompt" muon (muon directly originating from LFV Higgs boson decay)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, -1, leg2P4.pt(), leg2P4.eta(), leg2P4.phi(), leg2P4.mass(), 10)); // tau -> tau_h nu (3-prong) decay
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

  svFitAlgo.enableLogM(3.);
  //svFitAlgo.disableLogM(); 
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  svFitAlgo.setLikelihoodFileName("testClassicSVfitLFV.root");

  svFitAlgo.integrate(measuredEvent);
  bool isValidSolution = svFitAlgo.isValidSolution();

  double mass = svFitAlgo.getHistogramAdapter()->getMass();
  double massErr = svFitAlgo.getHistogramAdapter()->getMassErr();
 
  if ( isValidSolution ) {
    std::cout << "found valid solution: mass = " << mass << " +/- " << massErr << " (expected value = 126.12 +/- 16.9431)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass - 126.12) / 126.12) > 0.001) return 1;
  if (std::abs((massErr - 16.9431) / 16.9431) > 0.001) return 1;
    
  return 0;
}
