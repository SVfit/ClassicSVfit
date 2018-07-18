#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TMath.h"
#include "TVector2.h"
#include "TMatrixD.h"

#include "TF1.h"
#include "Math/BasicMinimizer.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::Likelihood(){

  covMET.ResizeTo(2,2);

  compnentsBitWord.reset();
  enableComponent(fastMTT::MASS);
  enableComponent(fastMTT::MET);
  enableComponent(fastMTT::PX);
  enableComponent(fastMTT::PY);
  //enableComponent(fastMTT::ENERGY);
  //enableComponent(fastMTT::IP);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::~Likelihood(){ }
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
 
  if(aLeg1DecayType==classic_svFit::MeasuredTauLepton::kTauToHadDecay && mVisLeg1>1.5){
    mVisLeg1 = 0.3;
  }  
  if(aLeg2DecayType==classic_svFit::MeasuredTauLepton::kTauToHadDecay && mVisLeg2>1.5){
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
void Likelihood::setIP3D(const double & aIP3DLeg1,
			 const double & aIP3DLeg2){

  ip3DLeg1 = aIP3DLeg1;
  ip3DLeg2 = aIP3DLeg2;
  
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
  if(parameters.size()<2) parameters = std::vector<double>{6, 1.0/1.15};

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::enableComponent(fastMTT::likelihoodComponent aCompIndex){

  compnentsBitWord.set(aCompIndex);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::disableComponent(fastMTT::likelihoodComponent aCompIndex){

  compnentsBitWord.reset(aCompIndex);
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
  if(delta<0){return std::make_tuple(0, 0); }

  double twoA = 2.0*(mVis2 + pVis2*sinGJ2);
  double solution1 = (b2 - sqrt(delta))/twoA;
  double solution2 = (b2 + sqrt(delta))/twoA;

  double tauEnergy1 = sqrt(pow(solution1,2) + mTau2);
  double tauEnergy2 = sqrt(pow(solution2,2) + mTau2);

  return std::make_tuple(tauEnergy1, tauEnergy2);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::energyLikelihood(const LorentzVector & tauP4,
				    const LorentzVector & visP4,
				    const double & cosGJ,
				    int decayMode) const{

  if(decayMode<10 || decayMode>14) return 1.0;
  
  std::tuple<double, double> energySolutions = energyFromCosGJ(visP4, cosGJ);
  if( std::get<0>(energySolutions)<1E-3 &&
      std::get<1>(energySolutions)<1E-3) return 1.0;

 
  std::vector<double> mean = {-0.05, 0.0};
  std::vector<double> sigma = {0.11, 0.08};

  //std::vector<double> mean = {-4.96169e-02, -4.96169e-02};
  //std::vector<double> sigma = {1.14987e-01,  1.14987e-01};

  ///Generator level cosGJ
  //std::vector<double> mean = {0.0, 0.0};
  //std::vector<double> sigma = {3.5e-02,  3.5e-02};
  
  double tauE = tauP4.E();
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

  double mpv = mean[1];  
  likelihoodSolution2 = TMath::Landau(pull2, mpv, sigma[1], true);
 
  double value = likelihoodSolution1 + likelihoodSolution2;
     
  return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::massLikelihood(const double & m) const{

  double coeff1 = parameters[0];
  double coeff2 = parameters[1];  
  double mScaled = m*coeff2;

  if(mScaled<mVis) return 0.0;

  const double & mTau = classic_svFit::tauLeptonMass;

  double x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  double x2Min = std::max(std::pow(mVisLeg2/mTau,2), std::pow(mVis/mScaled,2));
  double x2Max = std::min(1.0, std::pow(mVis/mScaled,2)/x1Min);
  if(x2Max<x2Min) return 0.0;

  double jacobiFactor = 2.0*std::pow(mVis,2)*std::pow(mScaled,-coeff1);
  double x2IntegralTerm = log(x2Max)-log(x2Min);
    
  double value = x2IntegralTerm;
  if(leg1DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    double mNuNuIntegralTermLeg1 = std::pow(mVis/mScaled,2)*(std::pow(x2Max,-1) - std::pow(x2Min,-1));
    value += mNuNuIntegralTermLeg1;
  }
  if(leg2DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    double mNuNuIntegralTermLeg2 = std::pow(mVis/mScaled,2)*x2IntegralTerm - (x2Max - x2Min);
    value += mNuNuIntegralTermLeg2;
  }
  
  ///The E9 factor to get values around 1.0
  value *=  1E8*jacobiFactor;
  
  return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::ptLikelihood(const double & pTTauTau, int type) const{

  ///Protection against numerical singularity in phase space volume.
  if(std::abs(pTTauTau)<0.5) return 0.0;

  const double & mTau = classic_svFit::tauLeptonMass;
  double pT1 = 0.0;
  double pT2 = 0.0;

  if(type==0){
    pT1 = leg1P4.Px();
    pT2 = leg2P4.Px();    
  }
  else if(type==1){
    pT1 = leg1P4.Py();
    pT2 = leg2P4.Py();
  }
  else{
    pT1 = leg1P4.Pz();
    pT2 = leg2P4.Pz();
  }

  Double_t x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  Double_t x2Min = std::min(1.0, std::pow(mVisLeg2/mTau,2));

   Double_t x2Max = 1.0;
   Double_t x1Max = 1.0;

   Double_t a_x2 = x1Min*pT2/(x1Min*pTTauTau - pT1);
   Double_t b_x2 = x1Max*pT2/(x1Max*pTTauTau - pT1);

   bool is_x1_vs_x2_falling = (-pT2*pT1)<0;
   bool x2_vs_x1_hasSingularity = pT1/pTTauTau>0.0 &&  pT1/pTTauTau<1.0;
   double x1_singularity = pT1/pTTauTau;
   if(x2_vs_x1_hasSingularity && x1_singularity<x1Min) return 0.0;


   if(is_x1_vs_x2_falling){
     x2Min = std::max(x2Min, b_x2);
     x2Max = std::min(x2Max, a_x2);
     if(x2_vs_x1_hasSingularity && x2Max<0) x2Max = 1.0;
   }
   else{
     x2Min = std::max(x2Min, a_x2);
     x2Max = std::min(x2Max, b_x2);
     if(x2_vs_x1_hasSingularity && x2Max<0) x2Max = 1.0;
   }

   if(x2Min<0) x2Min = 0.0;
   if(x2Min>x2Max) return 0.0;

  Double_t mNuNuIntegral = 0.0;
  Double_t x2 = std::min(1.0, x2Max);

  Double_t term1 = pT2-pTTauTau*x2;
  Double_t log_term1 = log(std::abs(term1));
  
  Double_t integralMax = (pT1*(pTTauTau*x2+pow(pT2,2)/term1+2*pT2*log_term1))/pow(pTTauTau,3);
  if(leg1DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    mNuNuIntegral =  -pow(pT1,2)*(2*pTTauTau*x2 + (pow(pT2,2)*(5*pT2 - 6*pTTauTau*x2))/pow(term1,2) +
				  6*pT2*log_term1)/(2*pow(pTTauTau,4));
  }
  if(leg2DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    mNuNuIntegral += -pT1/(2*pow(pTTauTau,5))*
                     (2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + 
		      pow(pTTauTau,2)*(-pT1 + pTTauTau)*pow(x2,2) +
		      (pow(pT2,4)*pT1)/pow(term1,2) + (2*pow(pT2,3)*(-4*pT1 + pTTauTau))/term1 + 
		      6*pow(pT2,2)*(-2*pT1 + pTTauTau)*log_term1); 
    }
  integralMax += mNuNuIntegral;
  
  x2 = x2Min;
  term1 = pT2-pTTauTau*x2;
  log_term1 = log(std::abs(term1));
  
  Double_t integralMin = (pT1*(pTTauTau*x2+pow(pT2,2)/term1+2*pT2*log_term1))/pow(pTTauTau,3);
  if(leg1DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    mNuNuIntegral =  -pow(pT1,2)*(2*pTTauTau*x2 + (pow(pT2,2)*(5*pT2 - 6*pTTauTau*x2))/pow(term1,2) +
				  6*pT2*log_term1)/(2*pow(pTTauTau,4));
  }
  if(leg2DecayType!=classic_svFit::MeasuredTauLepton::kTauToHadDecay){
  mNuNuIntegral += -pT1/(2*pow(pTTauTau,5))*
                   (2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + 
		    pow(pTTauTau,2)*(-pT1 + pTTauTau)*pow(x2,2) +
		    (pow(pT2,4)*pT1)/pow(term1,2) + (2*pow(pT2,3)*(-4*pT1 + pTTauTau))/term1 + 
		    6*pow(pT2,2)*(-2*pT1 + pTTauTau)*log_term1);  
  }

  integralMin += mNuNuIntegral;

  Double_t value  = integralMax - integralMin;

  if(std::abs(value)>100){
    std::cout<<"value: "<<value
	     <<" pTTauTau: "<<pTTauTau
	     <<std::endl;
  }
  
  ///The 1E4 factor to get values around 1.0
  value *= 1E4;

  return std::abs(value);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::ip3DLikelihood(const LorentzVector & tauP4,
				  const double & sinGJ,
				  const double & ip3D) const{

  const double & cTauLifetime = classic_svFit::cTauLifetime;
  const double & tauLeptonMass = classic_svFit::tauLeptonMass;
  
  double gammaBeta = tauP4.P()/tauLeptonMass;

  double arg = gammaBeta*cTauLifetime*std::abs(sinGJ);
  double value = 1.0/arg*exp(-ip3D/arg);

  return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
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

  const double & mTau = classic_svFit::tauLeptonMass;
  double x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  double x2Min = std::min(1.0, std::pow(mVisLeg2/mTau,2));
  
  if(x[0]<x1Min || x[1]<x2Min) return 0.0;
  
  testP4 = leg1P4*(1.0/x[0]) + leg2P4*(1.0/x[1]);
  testMET = testP4 - leg1P4 - leg2P4;

  double value = -1.0;
  if(compnentsBitWord.test(fastMTT::MET)) value *= metTF(recoMET, testMET, covMET);
  if(compnentsBitWord.test(fastMTT::MASS)) value *= massLikelihood(testP4.M());    
  if(compnentsBitWord.test(fastMTT::PX)) value *= ptLikelihood(testP4.Px(), 0);
  if(compnentsBitWord.test(fastMTT::PY)) value *= ptLikelihood(testP4.Py(), 1);
  if(compnentsBitWord.test(fastMTT::ENERGY)){
    value *= energyLikelihood((leg1P4*(1.0/x[0])),
			      leg1P4, cosGJLeg1, leg1DecayMode);

    value *= energyLikelihood((leg2P4*(1.0/x[1])),
			      leg2P4, cosGJLeg2, leg2DecayMode);
    
  }
  if(compnentsBitWord.test(fastMTT::IP)){    
    double sinGJ = sqrt(1.0 - cosGJLeg1*cosGJLeg1);
    value *= ip3DLikelihood(leg1P4*(1.0/x[0]),
			    sinGJ, ip3DLeg1);
    
    sinGJ = sqrt(1.0 - cosGJLeg2*cosGJLeg2);
    value *= ip3DLikelihood(leg2P4*(1.0/x[1]),
			    sinGJ, ip3DLeg2);
  }  
  return value;
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
  minimumPosition = initialValues;
  minimumValue = 999.0;

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
void FastMTT::enableComponent(fastMTT::likelihoodComponent aCompIndex){

  myLikelihood.enableComponent(aCompIndex);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::disableComponent(fastMTT::likelihoodComponent aCompIndex){

  myLikelihood.disableComponent(aCompIndex);

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
  myLikelihood.setIP3D(aLepton1.ip3D(), aLepton2.ip3D());
  myLikelihood.setMETInputs(aMET, covMET);

  scan();
  //minimize();

  tau1P4 = aLepton1.p4()*(1.0/minimumPosition[0]);
  tau2P4 = aLepton2.p4()*(1.0/minimumPosition[1]);
  bestP4 = tau1P4 + tau2P4;
  
  if(aLepton1.type() != classic_svFit::MeasuredTauLepton::kTauToHadDecay ||
     aLepton2.type() != classic_svFit::MeasuredTauLepton::kTauToHadDecay){
    bestP4 *= 1.03;
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

  minimizer->SetVariableLimits(0, 0.01, 1.0);
  minimizer->SetVariableLimits(1, 0.01, 1.0);

  minimizer->SetVariable(0, "x1", 0.5, 0.1);
  minimizer->SetVariable(1, "x2", 0.5, 0.1);
    
  minimizer->Minimize();

  const double *theMinimum = minimizer->X();
  minimumPosition[0] = theMinimum[0];
  minimumPosition[1] = theMinimum[1];
  minimumValue = minimizer->MinValue();

   if(minimizer->Status()!=0){
     std::cout<<" minimizer "
	      <<" Status: "<<minimizer->Status()
	      <<" nCalls: "<<minimizer->NCalls()
	      <<" nIterations: "<<minimizer->NIterations()
	      <<" x1Max: "<<theMinimum[0]
	      <<" x2Max: "<<theMinimum[1]
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
  double bestLH = 0.0;

  double x[2] = {0.5, 0.5};
  double theMinimum[2] = {0.75, 0.75};
  int nGridPoints = 100;
  int nCalls = 0;
  for(int iX2 = 1; iX2<nGridPoints;++iX2){
    x[1] = 1.0*(double)iX2/nGridPoints;
    for(int iX1 = 1; iX1<nGridPoints;++iX1){
      x[0] = 1.0*(double)iX1/nGridPoints;

      lh = myLikelihood.value(x);

      ++nCalls;
      if(lh<bestLH){
	bestLH = lh;
	theMinimum[0] = x[0];
	theMinimum[1] = x[1];
      }
    }
  }
  
  minimumPosition[0] = theMinimum[0];
  minimumPosition[1] = theMinimum[1];
  minimumValue = bestLH;

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

  return std::make_tuple(minimumPosition[0], minimumPosition[1]);
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getBestLikelihood() const{

  return minimumValue;
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getLikelihoodForX(double *x) const{

   return myLikelihood.value(x);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
