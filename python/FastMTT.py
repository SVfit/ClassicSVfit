import numpy as np
from vector import obj, Vector4D
from MeasuredTauLepton import *

###Reference: https://github.com/SVfit/ClassicSVfit/blob/fastMTT_2024/src/FastMTT.cc ###

class Likelihood:
    def __init__(self, enable_MET = True, enable_mass = True, enable_px = True, enable_py = True):
        #setMETinputs
        self.recoMET = Vector4D(px = 0.0, py = 0.0, pz = 0.0, E = 0.0)
        self.covMET = np.ones(2, 2)

        #setParameters
        self.coeff1 = 6
        self.coeff2 = 1.15

        #setLeptonInputs
        self.leg1P4 = Vector4D(px = 0.0, py = 0.0, pz = 0.0, E = 0.0)
        self.leg2P4 = Vector4D(px = 0.0, py = 0.0, pz = 0.0, E = 0.0)
        
        #visible invariant mass
        #eq. (4)
        self.mvis = 0.0

        self.mvisleg1 = 0.0
        self.mvisleg2 = 0.0

        self.mVisOverTauSquare1 = 0.0
        self.mVisOverTauSquare2 = 0.0
            
        self.mTau = 0.0
        self.leg1DecayType = 0
        self.leg2DecayType = 0
        self.leg1DecayMode = 0
        self.leg2DecayMode = 0

        self.enable_MET = enable_MET
        self.enable_mass = enable_mass
        self.enable_px = enable_px
        self.enable_py = enable_py

        return
    

    def setMETinputs(self, aMET, aCovMET):
        self.recoMET = aMET
        self.covMET = aCovMET

    def setParameters(self, aPars):
        self.coeff1 = aPars[0]
        self.coeff2 = aPars[1]

    def setLeptonInputs(self, aLeg1P4, aLeg2P4, aLeg1DecayType, aLeg2DecayType, aLeg1DecayMode, aLeg2DecayMode, tauLeptonMass = 1.777):
        self.leg1P4 = aLeg1P4
        self.leg2P4 = aLeg2P4
        
        #visible invariant mass
        #eq. (4)
        self.mvis = (self.leg1P4 + self.leg2P4).mass()

        self.mvisleg1 = self.leg1P4.mass()
        self.mvisleg2 = self.leg2P4.mass()

        self.mVisOverTauSquare1 = (self.mvisleg1/self.mTau)**2
        self.mVisOverTauSquare2 = (self.mvisleg2/self.mTau)**2

        if aLeg1DecayType==MeasuredTauLepton.kTauToHadDecay and self.mvis1>1.5:
            self.mvisleg1 = 0.3
        if aLeg2DecayType==MeasuredTauLepton.kTauToHadDecay and self.mvis2>1.5:
            self.mvisleg2 = 0.3
            
        self.mTau = tauLeptonMass
        self.leg1DecayType = aLeg1DecayType
        self.leg2DecayType = aLeg2DecayType
        self.leg1DecayMode = aLeg1DecayMode
        self.leg2DecayMode = aLeg2DecayMode


    def massLikelihood(self, m: float):
        mScaled = m*self.coeff2
        if mScaled<self.mvis:
            return 0.0
        mVS2 = (self.mvis/mScaled)**2
        x1Min = min(1.0, )
        x2Min = max(self.mVisOverTauSquare2, mVS2)
        x2Max = min(1.0, mVS2/x1Min)
        
        if x2Min>x2Max:
            return 0.0
        
        jacobiFactor = 2.0*self.mvis**2*mScaled**(-self.coeff1)
        x2IntegralTerm = np.log(x2Max/x2Min)

        value = x2IntegralTerm

        if self.leg1DecayType != MeasuredTauLepton.kTauToHadDecay:
            value += mVS2*x2IntegralTerm - (x2Max - x2Min)

        value *= 1E9*jacobiFactor
        return value
    
    def ptLikelihood(self, pTTauTau: float, type: int):
        if np.abs(pTTauTau)<0.5:
            return 0.0
        
        if type == 0:
            pT1 = self.leg1P4.px()
            pT2 = self.leg2P4.px()
        elif type == 1:
            pT1 = self.leg1P4.py()
            pT2 = self.leg2P4.py()
        elif type == 2:
            pT1 = self.leg1P4.pz()
            pT2 = self.leg2P4.pz()

        x1Min = np.min(1.0, self.mVisOverTauSquare1)
        x2Min = np.min(1.0, self.mVisOverTauSquare2)

        x1Max = 1.0
        x2Max = 1.0

        a_x2 = x1Min *pT2/(x1Min*pTTauTau - pT1)
        b_x2 = x1Max*pT2/(x1Max*pTTauTau - pT1)

        x1_singularity = pT1/pTTauTau
        x2_vs_x1_singularity = x1_singularity>0.0 and x1_singularity<1.0
        if x2_vs_x1_singularity and x1_singularity<x1Min:
            return 0.0

        if (-pT2*pT1)<0:
            x2Min = np.max(x2Min, b_x2)
            x2Max = np.min(x2Max, a_x2)
            if x2_vs_x1_singularity and x2Max<0:
                x2Max = 1.0
        else:
            x2Min = np.max(x2Min, a_x2)
            x2Max = np.min(x2Max, b_x2)
            if x2_vs_x1_singularity and x2Max<0:
                x2Max = 1.0

        if x2Min<0:
            x2Min = 0.0
        
        if x2Min > x2Max:
            return 0.0
        
        mNuNuIntegral = 0.0
        x2 = np.min(1.0, x2Max)

        term1 = pT2 - pTTauTau*x2
        log_term1 = np.log(np.abs(term1))

        integralMax = pT1*(pTTauTau*x2 + pT2**2/term1 + 2*pT2*log_term1)/pTTauTau**3

        if self.leg1DecayType != MeasuredTauLepton.kTauToHadDecay:
            mNuNuIntegral = -pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term1**2 + 6*pT2*log_term1)/(2*pTTauTau**4)

        if self.leg2DecayType != MeasuredTauLepton.kTauToHadDecay:
            mNuNuIntegral += -pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term1**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term1 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term1)

        integralMax += mNuNuIntegral

        x2 = x2Min
        term2 = pT2 - pTTauTau*x2
        log_term2 = np.log(np.abs(term2))

        integralMin = pT1*(pTTauTau*x2+pT2**2/term2+2*pT2*log_term2)/pTTauTau**3
        if self.leg1DecayType != MeasuredTauLepton.kTauToHadDecay:
            mNuNuIntegral = -pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term2**2+6*pT2*log_term2)/(2*pTTauTau**4)

        if self.leg2DecayType != MeasuredTauLepton.kTauToHadDecay:
            mNuNuIntegral += -pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term2**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term2 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term2)

        integralMin += mNuNuIntegral

        value = integralMax - integralMin
        value*=1E4

        return np.abs(value)
    
    def metTF(metP4: Vector4D, nuP4: Vector4D, covMET):
        aMETx = metP4.x
        aMETy = metP4.y

        covDET = np.det(covMET)
        
        if covDET < 1E-10:
            print(f"Error: Cannot invert MET covariance matrix (det=0)! aMETx: {aMETx}, aMETy: {aMETy}")
            return 0.0

        constMET = 1/2/np.pi/np.sqrt(covDET)
        residualX = aMETx - nuP4.x
        residualY = aMETy - nuP4.y

        pull2 = residualX*(covMET[0][0]*residualX + covMET[0][1]*residualY) + residualY*(covMET[1][0]*residualX + covMET[1][1]*residualY)

        return constMET*np.exp(-0.5*pull2)
    
    #Uwaga#
    #ComponentParams zastąpią enable/disable Component
    #które (jeśli dobrze rozumiem) uwzględniają lub nie konkretne prawdopodobieństwa (pt, metTF, mass) w końcowym wyniku
    #Możemy to zrobić w value albo w __init__ - jeszcze nie wiem, co jest potrzebne do poprawnego działania kodu

    def value(self, x, ComponentsParams = True):
        x1Min = np.min(1.0, self.mVisOverTauSquare1)
        x2Min = np.min(1.0, self.mVisOverTauSquare2)

        if x[0]<x1Min or x[1]<x2Min:
            return 0.0
        
        testP4 = self.leg1P4/x[0]+self.leg2P4/x[1]
        testMET = testP4 - self.leg1P4 - self.leg2P4

        value = -1.0
        if self.enable_MET:
            value *= self.metTF(self.recoMET, testMET, self.covMET)
        if self.enable_mass:
            value *= self.massLikelihood(testP4.mass())
        if self.enable_px:
            value *= self.ptLikelihood(testP4.pt(), 0)
        if self.enable_py:
            value *= self.ptLikelihood(testP4.pt(), 1)

        return value


        #Implementacja
        #zawiera odwołania do massLikelihood, ptLikelihood, metTF
        #o ile są odblokowane przy definicji funkcji
        return

###UWAGA###
#Ponieważ da się to jednak prosto napisać bez dziedziczenia funkcji, to spróbujemy zarówno z numbą, jak i innymi plikami jit

class FastMTT(Likelihood):
    def __init__(self):
        self.myLikelihood = Likelihood()
        self.BestLikelihood = 0.0
        self.BestX = np.array([0.0, 0.0])
        self.bestP4 = 0.0
        #inicjalizacja
        #parametry aPars i ComponentParams przez dziedziczenie z Likelihood
        #Odpuścimy inicjalizację rzeczy do minimalizacji, skoro i tak jej nie ma w kodzie
        return
    
    def run(self, measuredTauLeptons: np.ndarray, measuredMETx, measuredMETy, covMET) -> np.ndarray:
        if measuredTauLeptons != 2:
            print(f"Number of MeasuredTauLepton is {measuredTauLeptons.size()}. A user shouls pass exactly two leptons.\n")
            return
        
        #Waiting for sorting#
        sortedMeasuredTauLeptons = measuredTauLeptons
        metLenght = np.sqrt(measuredMETx**2 + measuredMETy**2)
        aMET = Vector4D(px = measuredMETx, py = measuredMETy, pz = 0.0, E = metLenght)
        aLepton1 = MeasuredTauLepton[0]
        aLepton2 = MeasuredTauLepton[1]

        self.myLikelihood.setMETinputs(aMET, covMET)
        self.myLikelihood.setLeptonInputs(aLepton1.p4(), aLepton2.p4(), aLepton1.type, aLepton2.type, aLepton1.decayMode, aLepton2.decayMode)

        self.scan()
        tau1P4 = aLepton1.p4()*(1/self.minimumPosition[0])
        tau2P4 = aLepton2.p4()*(1/self.minimumPosition[1])
        self.bestP4 = tau1P4 + tau2P4
    
    def compareLeptons(self, measuredTauLepton1: MeasuredTauLepton, measuredTauLepton2: MeasuredTauLepton): #używane w run
        if (measuredTauLepton1.type == MeasuredTauLepton.kTauToElecDecay or measuredTauLepton1.type == MeasuredTauLepton.kTauToMuDecay) and measuredTauLepton2.type == MeasuredTauLepton.kTauToHadDecay:
            return True
        if (measuredTauLepton2.type == MeasuredTauLepton.kTauToElecDecay or measuredTauLepton2.type == MeasuredTauLepton.kTauToMuDecay) and measuredTauLepton1.type == MeasuredTauLepton.kTauToHadDecay:
            return False
        return measuredTauLepton1.pt > measuredTauLepton2.pt
    
    def myLikelihoodValue(self, x):
        return self.myLikelihood.value(x)
    
    def scan(self): #używane w run
        lh = 0.0
        bestLH = 0.0
        x = np.array([0.5, 0.5])
        theMinimum = np.array([0.75, 0.75])
        nGridPoints = 100
        gridFactor = 1.0/nGridPoints
        nCalls = 0
        for iX2 in range(1, nGridPoints):
            x[1] = iX2*gridFactor
            for iX1 in range(1, nGridPoints):
                x[0] = iX1*gridFactor
                lh = self.myLikelihood.value(x)
                nCalls += 1
                if lh < bestLH:
                    bestLH = lh
                    theMinimum[0] = x[0]
                    theMinimum[1] = x[1]

        self.minimumPosition[0] = theMinimum[0]
        self.minimumPosition[1] = theMinimum[1]
        self.minimumValue = bestLH

        #Używa MyLikelihood.value
        #Implementacja z pętlami
        #lub gradient_descent (w artykule było wspomniane, że może poprawić szybkość algorytmu)
        return
    
#Do wszystkiego dołożymy do testowania funkcje do pomiaru czasu (np. import time), tak jak w oryginalnym kodzie