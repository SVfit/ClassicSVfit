import numpy as np
from vector import obj, Vector4D
from MeasuredTauLepton import *

###Reference: https://github.com/SVfit/ClassicSVfit/blob/fastMTT_2024/src/FastMTT.cc ###

class Likelihood:
    def __init__(self, aLeg1P4: Vector4D, aLeg2P4: Vector4D, aLeg1DecayType, aLeg2DecayType, aLeg1DecayMode, aLeg2DecayMode, aMET, aCovMET, aPars = np.array([6, 1/1.15]), enable_MET = True, enable_mass = True, enable_px = True, enable_py = True):
        #setMETinputs
        self.recoMET = aMET
        self.covMET = aCovMET

        #setParameters
        self.coeff1 = aPars[0]
        self.coeff2 = aPars[1]

        #setLeptonInputs
        self.leg1P4 = aLeg1P4
        self.leg2P4 = aLeg2P4
        sum_of_legs = self.leg1P4 + self.leg2P4
        
        #visible invariant mass
        #eq. (4)
        self.mvis = (self.leg1P4 + self.leg2P4).mass()

        self.mvisleg1 = self.leg1P4.mass()
        self.mvisleg2 = self.leg2P4.mass()

        self.mVisOverTauSquare1 = (self.mvisleg1/self.mTau)**2
        self.mVisOverTauSquare2 = (self.mvisleg2/self.mTau)**2

        if aLeg1DecayType==MeasuredTauLepton.kTauToHadDecay and self.mvis1>1.5:
            self.mvis1 = 0.3
        if aLeg2DecayType==MeasuredTauLepton.kTauToHadDecay and self.mvis2>1.5:
            self.mvis2 = 0.3
            
        self.mTau = tauLeptonMass
        self.leg1DecayType = aLeg1DecayType
        self.leg2DecayType = aLeg2DecayType
        self.leg1DecayMode = aLeg1DecayMode
        self.leg2DecayMode = aLeg2DecayMode

        self.enable_MET = enable_MET
        self.enable_mass = enable_mass
        self.enable_px = enable_px
        self.enable_py = enable_py

        return
    
    def massLikelihood(self, m):
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
    
    def ptLikelihood(self, pTTauTau, type):
        if np.abs(pTTauTau)<0.5:
            return 0.0
        x1Min = np.min(1.0, self.mVisOverTauSquare1)
        x2Min = np.min(1.0, self.mVisOverTauSquare2)

        x1Max = 1.0
        x2Max = 1.0

        #Uzupełnić z PowerTable
        a_x2 = 0.0
        a_x1 = 0.0
        


        #Implementacja
        return #value
    
    def metTF(metP4: Vector4D, nuP4: Vector4D, covMET):
        aMETx = metP4.px
        aMETy = metP4.py

        covDET = np.det(covMET)
        
        if covDET < 1E-10:
            print(f"Error: Cannot invert MET covariance matrix (det=0)! aMETx: {aMETx}, aMETy: {aMETy}")
            return 0.0

        constMET = 1/2/np.pi/np.sqrt(covDET)
        residualX = aMETx - nuP4.px
        residualY = aMETy - nuP4.py

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
#numba nie obsługuje dziedziczenia funkcji
#więc preferowanym kompilatorem będzie chyba jednak jit
#Tym niemniej przepisanie kodu do jednej klasy też nie powinno być ewentualnie problemem (prędzej dla MeasuredTauLepton niż dla likelihood)
#jeśli performance z numbą byłby znacznie większy
#więc domyślnie napisałbym kod z jit, a potem ew. sprawdził jak to działa z numbą
#(jeśli uznamy, że implementacja nie jest zbyt ciężka w porównaniu do wyniku)

class FastMTT(Likelihood):
    def __init__(self):
        #inicjalizacja
        #parametry aPars i ComponentParams przez dziedziczenie z Likelihood
        #Odpuścimy inicjalizację rzeczy do minimalizacji, skoro i tak jej nie ma w kodzie
        return
    
    def run(self, measuredTauLeptons: np.ndarray, measuredMETx, measuredMETy, covMET) -> np.ndarray:
        return
    
    def compareLeptons(self, measuredTauLepton1: MeasuredTauLepton, measuredTauLepton2: MeasuredTauLepton): #używane w run
        #Implementacja
        return
    
    def scan(self): #używane w run
        #Używa MyLikelihood.value
        #Implementacja z pętlami
        #lub gradient_descent (w artykule było wspomniane, że może poprawić szybkość algorytmu)
        return
    
#Do wszystkiego dołożymy do testowania funkcje do pomiaru czasu (np. import time), tak jak w oryginalnym kodzie