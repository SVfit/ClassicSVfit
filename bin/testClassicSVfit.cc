
/**
   \class testClassicSVfit testClassicSVfit.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm
*/

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"

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
	//Run FastMTT
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
	LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	std::cout<<std::endl;
	std::cout << "FastMTT found best p4 with mass = " << ttP4.M()
			<< " (expected value = 108.991),"
			<<std::endl;
	std::cout<<"Real Time =   "<<aFastMTTAlgo.getRealTime("scan")<<" seconds "
		 <<" Cpu Time =   "<<aFastMTTAlgo.getCpuTime("scan")<<" seconds"<<std::endl;
	if(std::abs(ttP4.M() -  108.991)>1E-3*108.991){
		std::cout<<"error value: "<<std::abs(ttP4.M() -  108.991)<<std::endl;
		return 1;
	}

  return 0;
}
