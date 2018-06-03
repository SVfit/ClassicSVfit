#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TMath.h"
#include "TVector2.h"
#include "TMatrixD.h"

#include "TF1.h"
#include "Math/BasicMinimizer.h"

#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::Likelihood(){

  covMET.ResizeTo(2,2);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::~Likelihood(){}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setLeptonInputs(const LorentzVector & aLeg1P4,
                                 const LorentzVector & aLeg2P4,
                                 int aLeg1DecayType, int aLeg2DecayType,
				 int aLeg1DecayMode, int aLeg2DecayMode){
  leg1P4 = aLeg1P4;
  leg2P4 = aLeg2P4;

  mVis = (leg1P4 + leg2P4).M();
  mVisLeg1 = leg1P4.M();
  mVisLeg2 = leg2P4.M();
  ///Setting fixed hadronic tau mass avoids
  ///cases when solution is not found.
  if(aLeg2DecayType==classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    mVisLeg2 = 0.3;
  }

  leg1DecayType = aLeg1DecayType;
  leg2DecayType = aLeg2DecayType;

  leg1DecayMode = aLeg1DecayMode;
  leg2DecayMode = aLeg2DecayMode;
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setCosGJ(const double & aCosGJLeg1,
			  const double & aCosGJLeg2){

  cosGJLeg1 = aCosGJLeg1;
  cosGJLeg2 = aCosGJLeg2;
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setMETInputs(const LorentzVector & aMET,
                              const TMatrixD& aCovMET){
  recoMET = aMET;
  covMET = aCovMET;

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setParameters(const std::vector<double> & aPars){

  parameters = aPars;

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
std::tuple<double, double> Likelihood::energyFromCosGJ(const LorentzVector & visP4,
						       const double & cosGJ) const{

  const double & mTau2 = classic_svFit::tauLeptonMass2;

  double mVis =  visP4.M();
  double mVis2 = std::pow(mVis,2);
  double pVis =  visP4.P();
  double pVis2 = std::pow(pVis,2);
  double cosGJ2 = std::pow(cosGJ,2);
  double sinGJ2 = 1.0 - cosGJ2;

  double b2 = (mVis2 + mTau2)*pVis*cosGJ;

  double delta = (mVis2 + pVis2)*(std::pow(mVis2 - mTau2, 2) - 4.0*mTau2*pVis2*sinGJ2);
  if(delta<0){
    sinGJ2 = 0;
    delta = (mVis2 + pVis2)*(std::pow(mVis2 - mTau2, 2) - 4.0*mTau2*pVis2*sinGJ2);
  }

  double twoA = 2.0*(mVis2 + pVis2*sinGJ2);
  double solution1 = (b2 - sqrt(delta))/twoA;
  double solution2 = (b2 + sqrt(delta))/twoA;

  double tauEnergy1 = sqrt(pow(solution1,2) + mTau2);
  double tauEnergy2 = sqrt(pow(solution2,2) + mTau2);

  return std::make_tuple(tauEnergy1, tauEnergy2);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::energyLikelihood(const double & tauE,
				    const LorentzVector & visP4,
				    const double & cosGJ,
				    int decayMode) const{

  if(decayMode<10 || decayMode>14) return 1.0;
  
  std::tuple<double, double> energySolutions = energyFromCosGJ(visP4, cosGJ);

  std::vector<double> mean = {-0.163, 0.114};
  std::vector<double> sigma = {0.146,  0.16};

  double pull1 = std::get<0>(energySolutions) - tauE;
  pull1 /= tauE;
  pull1 -= mean[0];

  double gaussNorm = 1.0/sigma[0]/sqrt(2.0*M_PI);
  double likelihoodSolution1 = gaussNorm*TMath::Exp(-0.5*std::pow(pull1/sigma[0], 2));

  double pull2 = std::get<1>(energySolutions) - tauE;
  pull2 /= tauE;
  pull2 -= mean[1];
  
  gaussNorm = 1.0/sigma[1]/sqrt(2.0*M_PI);
  double likelihoodSolution2 = gaussNorm*TMath::Exp(-0.5*std::pow(pull2/sigma[1], 2));
  /*
  std::cout<<"tauE: "<<tauE
           <<" solution1: "<<std::get<0>(energySolutions)
	   <<" solution2: "<<std::get<1>(energySolutions)
	   <<" pull1: "<<pull1
    	   <<" pull2: "<<pull2
	   <<" likelihoodSolution1: "<<likelihoodSolution1
    	   <<" likelihoodSolution2: "<<likelihoodSolution2
	   <<std::endl;
  */

  if(tauE>std::get<0>(energySolutions) || tauE<0.8*std::get<0>(energySolutions)) return 0.0;
  else return likelihoodSolution1;
 
  return likelihoodSolution1;
  return std::max(likelihoodSolution1, likelihoodSolution2);
  
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::massLikelihood(const double & m) const{

  double coeff1 = parameters[0];
  double coeff2 = parameters[1];
  double mShift = m*coeff2;

  if(mShift<mVis) return 0.0;

  double mTau = 1.77685;

  double x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  double x2Min = std::max(std::pow(mVisLeg2/mTau,2), std::pow(mVis/mShift,2));
  double x2Max = std::min(1.0, std::pow(mVis/mShift,2)/x1Min);

  double jacobiFactor = 2.0*std::pow(mVis,2)*std::pow(mShift,-coeff1);
  double x2IntegralTerm = log(x2Max)-log(x2Min);
    
  double value = x2IntegralTerm;
  if(leg1DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    double mNuNuIntegralTermLeg1 = std::pow(mVis/mShift,2)*(std::pow(x2Max,-1) - std::pow(x2Min,-1));
    value += mNuNuIntegralTermLeg1;
  }
  if(leg2DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    double mNuNuIntegralTermLeg2 = std::pow(mVis/mShift,2)*x2IntegralTerm - (x2Max - x2Min);
    value += mNuNuIntegralTermLeg2;
  }
    
  ///The E12 factor to get values around 1.0
  value *=  1E12*jacobiFactor;

  return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::ptLikelihood(const double & pTTauTau, int type) const{

  double mTau = 1.77685;

  double pT1 = 0.0;
  double pT2 = 0.0;

  if(type==0){
    pT1 = leg1P4.Px();
    pT2 = leg2P4.Px();
  }
  else{
    pT1 = leg1P4.Py();
    pT2 = leg2P4.Py();
  }
   
  Double_t x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  Double_t x1Max = pT1/(pTTauTau - pT2);
  if(x1Max<0.0) x1Max = 1.0;
  if(x1Min>x1Max) return 0.0;
  
  Double_t x1 = std::min(1.0, x1Max);
  Double_t integralMax = (pT2*(pTTauTau*x1+pow(pT1,2)/(pT1-pTTauTau*x1)+2*pT1*log(std::abs(pT1-pTTauTau*x1))))/pow(pTTauTau,3);
  
  x1 = x1Min;
  Double_t integralMin = (pT2*(pTTauTau*x1+pow(pT1,2)/(pT1-pTTauTau*x1)+2*pT1*log(std::abs(pT1-pTTauTau*x1))))/pow(pTTauTau,3);

  Double_t value  = integralMax - integralMin;
  
  ///The E4 factor to get values around 1.0
  value *= 1E4;
  
  return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////x
double Likelihood::metTF(const LorentzVector & metP4,
                         const LorentzVector & nuP4,
                         const TMatrixD& covMET) const{

  double  aMETx = metP4.X();
  double  aMETy = metP4.Y();

  double invCovMETxx = covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy = covMET(0,0);
  double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet)<1E-10){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
	      <<"METx: "<<aMETy<<" METy: "<<aMETy
	      << std::endl;
    return 0;
  }
  double const_MET = 1./(2.*M_PI*TMath::Sqrt(covDet));

  double residualX = aMETx - (nuP4.X());
  double residualY = aMETy - (nuP4.Y());

  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
    residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2/=covDet;

  return const_MET*TMath::Exp(-0.5*pull2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double Likelihood::value(const double *x) const{

  LorentzVector testP4 = leg1P4*(1.0/x[0]) + leg2P4*(1.0/x[1]);
  LorentzVector testMET = testP4 - leg1P4 - leg2P4;
 
  double value = metTF(recoMET, testMET, covMET);
  //value *= massLikelihood(testP4.M());
  value *= ptLikelihood(testP4.Px(), 0);//TEST
  value *= ptLikelihood(testP4.Py(), 1);//TEST
  
  double llh1 = energyLikelihood((leg1P4*(1.0/x[0])).E(),
				 leg1P4, cosGJLeg1, leg1DecayMode);

  double llh2 = energyLikelihood((leg2P4*(1.0/x[1])).E(),
				 leg2P4, cosGJLeg2, leg2DecayMode);
  
  return -value;
  //return -value*llh1*llh2;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
FastMTT::FastMTT(){

  minimizerName = "Minuit2";
  minimizerAlgorithm = "Migrad";
  initialize();
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
FastMTT::~FastMTT(){

  delete minimizer;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::initialize(){

  minimizer = ROOT::Math::Factory::CreateMinimizer(minimizerName, minimizerAlgorithm);
  minimizer->SetMaxFunctionCalls(100000);
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(0.01);

  std::vector<std::string> varNames = {"x1", "x2"};
  nVariables = varNames.size();
  std::vector<double> initialValues(nVariables,0.5);
  std::vector<double> stepSizes(nVariables, 0.01);
  minimalizationResult = initialValues;

  for(unsigned int iVar=0; iVar<nVariables; ++iVar){
    minimizer->SetVariable(iVar, varNames[iVar].c_str(), initialValues[iVar], stepSizes[iVar]);
  }

  std::vector<double> shapeParams = {6, 1.0/1.15};
  setLikelihoodParams(shapeParams);
  likelihoodFunctor = new ROOT::Math::Functor(&myLikelihood, &Likelihood::value, nVariables);
  minimizer->SetFunction(*likelihoodFunctor);

  verbosity = 0;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::setLikelihoodParams(const std::vector<double> & aPars){

   myLikelihood.setParameters(aPars);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
bool FastMTT::compareLeptons(const classic_svFit::MeasuredTauLepton& measuredTauLepton1,
			     const classic_svFit::MeasuredTauLepton& measuredTauLepton2){

  using namespace classic_svFit;
  
  if ( (measuredTauLepton1.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton1.type() == MeasuredTauLepton::kTauToMuDecay) &&
       measuredTauLepton2.type() == MeasuredTauLepton::kTauToHadDecay  ) return true;
  if ( (measuredTauLepton2.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2.type() == MeasuredTauLepton::kTauToMuDecay) &&
       measuredTauLepton1.type() == MeasuredTauLepton::kTauToHadDecay ) return false;
  return ( measuredTauLepton1.pt() > measuredTauLepton2.pt() );
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::run(const std::vector<classic_svFit::MeasuredTauLepton>& measuredTauLeptons,
		  const double & measuredMETx, const double & measuredMETy,
		  const TMatrixD& covMET){

  bestP4 = LorentzVector();  
  ////////////////////////////////////////////
  
  if(measuredTauLeptons.size()!=2){
    std::cout<<"Number of MeasuredTauLepton is "<<measuredTauLeptons.size()
	     <<" a user shouls pass exactly two leptons."<<std::endl;
    return;
  }

  std::vector<classic_svFit::MeasuredTauLepton> sortedMeasuredTauLeptons = measuredTauLeptons; 
  std::sort(sortedMeasuredTauLeptons.begin(),
  	    sortedMeasuredTauLeptons.end(),
  	    compareLeptons);
 
  double metLength = sqrt(std::pow(measuredMETx, 2) +
			  std::pow(measuredMETy, 2));
  LorentzVector aMET =  LorentzVector(measuredMETx, measuredMETy, 0, metLength);


  const classic_svFit::MeasuredTauLepton & aLepton1 = measuredTauLeptons[0];
  const classic_svFit::MeasuredTauLepton & aLepton2 = measuredTauLeptons[1];
  
  myLikelihood.setLeptonInputs(aLepton1.p4(), aLepton2.p4(),
			       aLepton1.type(), aLepton2.type(),
			       aLepton1.decayMode(), aLepton2.decayMode());
  
  myLikelihood.setCosGJ(aLepton1.cosGJ(), aLepton2.cosGJ());
  myLikelihood.setMETInputs(aMET, covMET);

  scan();

  bestP4 = aLepton1.p4()*(1.0/minimalizationResult[0]) +
           aLepton2.p4()*(1.0/minimalizationResult[1]);

  if(aLepton1.type() != classic_svFit::MeasuredTauLepton::kTauToHadDecay &&
     aLepton2.type() != classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    bestP4 *= 1.036;
  }
  if(aLepton1.type() == classic_svFit::MeasuredTauLepton::kTauToHadDecay &&
     aLepton2.type() == classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    bestP4 *= 0.98;
  }  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::minimize(){

  clock.Reset();
  clock.Start("minimize");

  minimizer->SetVariableLimits(0, 0.0 ,1.0);
  minimizer->SetVariableLimits(1, 0.0 ,1.0);
  minimizer->Minimize();

  const double *theMinimum = minimizer->X();
  minimalizationResult[0] = theMinimum[0];
  minimalizationResult[1] = theMinimum[1];

   if(minimizer->Status()!=0){
     std::cout<<" minimizer "
	      <<" Status: "<<minimizer->Status()
	      <<" nCalls: "<<minimizer->NCalls()
	      <<" nIterations: "<<minimizer->NIterations()
	      <<" x1Max: "<<theMinimum[0]
	      <<" x2Max: "<<theMinimum[1]
	      <<" x3Max: "<<theMinimum[2]
	      <<" max LLH: "<<minimizer->MinValue()
	      <<" m: "<<bestP4.M()
	      <<std::endl;
}
  clock.Stop("minimize");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::scan(){

  clock.Reset();
  clock.Start("scan");

  double lh = 0.0;
  double maxLH = -99.0;

  double x[2] = {0.5, 0.5};
  double theMinimum[2] = {0.5, 0.5};  
  int nGridPoints = 100;
  int nCalls = 0;
  for(int iX2 = 1; iX2<nGridPoints;++iX2){
    x[1] = (double)iX2/nGridPoints;
    for(int iX1 = 1; iX1<nGridPoints;++iX1){
      x[0] = (double)iX1/nGridPoints;

      lh = - myLikelihood.value(x);
      ++nCalls;

      if(lh>maxLH){
	maxLH = lh;
	theMinimum[0] = x[0];
	theMinimum[1] = x[1];
      }
    }
  }
  minimalizationResult[0] = theMinimum[0];
  minimalizationResult[1] = theMinimum[1];

  clock.Stop("scan");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getCpuTime(const std::string & method){

  return clock.GetCpuTime(method.c_str());
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getRealTime(const std::string & method){

  return clock.GetRealTime(method.c_str());
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
std::tuple<double, double> FastMTT::getBestX() const{

  return std::make_tuple(minimalizationResult[0], minimalizationResult[1]);
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
