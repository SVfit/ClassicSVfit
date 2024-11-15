import numpy as np
import matplotlib.pyplot as plt
from MeasuredTauLepton import *
import time

#import jax
#import jax.numpy as jnp
#from jax import grad
#from jax.scipy.optimize import minimize as jax_minimize

#from numba import jit
#from scipy.optimize import minimize

###Main reference: https://github.com/SVfit/ClassicSVfit/blob/fastMTT_2024/src/FastMTT.cc ###


#Invariant mass calculation
def InvariantMass(aP4):
    energy_squared = aP4[..., 3]**2
    momentum_squared = aP4[..., 0]**2 + aP4[..., 1]**2 + aP4[..., 2]**2
    return np.sqrt(energy_squared - momentum_squared)

class Likelihood:
    def __init__(self, enable_MET = True, enable_mass = True, enable_px = False, enable_py = False):
        #METinputs
        self.recoMET = np.array([0.0, 0.0, 0.0, 0.0])
        self.covMET = np.ones((2, 2))

        #setParameters
        self.coeff1 = 6
        self.coeff2 = 1/1.15

        #LeptonInputs
        self.leg1P4 = np.array([0.0, 0.0, 0.0, 0.0])
        self.leg2P4 = np.array([0.0, 0.0, 0.0, 0.0])

        #Visible mass of both leptons
        self.mvis = np.array([0.0, 0.0, 0.0, 0.0])

        #Invariant mass of each lepton
        self.mvisleg1 = np.array([0.0])
        self.mvisleg2 = np.array([0.0])

        self.mVisOverTauSquare1 = np.array([0.0])
        self.mVisOverTauSquare2 = np.array([0.0])
         
        self.mTau = tauLeptonMass
        
        self.leg1DecayType = np.array([0.0])
        self.leg2DecayType = np.array([0.0])
        self.leg1DecayMode = np.array([0.0])
        self.leg2DecayMode = np.array([0.0])


        #Enable/disable likelihood channel
        self.enable_MET = enable_MET
        self.enable_mass = enable_mass

        #These are experimental and not used by main code
        self.enable_px = enable_px
        self.enable_py = enable_py

        return

    def setParameters(self, aPars):
        self.coeff1 = aPars[0]
        self.coeff2 = aPars[1]

    def setLeptonInputs(self, aLeg1P4, aLeg2P4, aLeg1DecayType, aLeg2DecayType, aLeg1DecayMode, aLeg2DecayMode):
        
        self.leg1DecayType = aLeg1DecayType
        self.leg2DecayType = aLeg2DecayType
        self.leg1DecayMode = aLeg1DecayMode
        self.leg2DecayMode = aLeg2DecayMode

        self.leg1P4 = aLeg1P4
        self.leg2P4 = aLeg2P4
        
        #visible invariant mass
        #eq. (4)
        self.mvis = InvariantMass(self.leg1P4 + self.leg2P4)
        
        self.mvisleg1[(aLeg1DecayType==1) & (self.mvisleg1>1.5)] = 0.3
        self.mvisleg2[(aLeg2DecayType==1) & (self.mvisleg2>1.5)] = 0.3

        self.mVisOverTauSquare1 = (self.mvisleg1/self.mTau)**2
        self.mVisOverTauSquare2 = (self.mvisleg2/self.mTau)**2


    def massLikelihood(self, m: np.ndarray):
        mScaled = m*self.coeff2

        mask1 = (mScaled < self.mvis[:, np.newaxis])
        
        mVS2 = (self.mvis[:, np.newaxis]/mScaled)**2
        
        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)
        x2Min = np.maximum(self.mVisOverTauSquare2[:, np.newaxis], mVS2)
        x2Max = np.minimum(1.0, mVS2/x1Min[:, np.newaxis])
        
        mask2 = (x2Min > x2Max)
        
        jacobiFactor = 2.0*self.mvis[:, np.newaxis]**2*mScaled**(-self.coeff1)
        x2IntegralTerm = np.log(x2Max/x2Min)

        value = x2IntegralTerm

        HadDecay1 = np.broadcast_to((self.leg1DecayType != 1)[:, np.newaxis], value.shape)
        value += HadDecay1 * mVS2 * (1 / x2Max - 1 / x2Min)

        HadDecay2 = np.broadcast_to((self.leg2DecayType != 1)[:, np.newaxis], value.shape)
        value += HadDecay2 * (mVS2*x2IntegralTerm - (x2Max - x2Min))

        value[mask1 | mask2] = 0.0

        value *= 1E9*jacobiFactor
        return value
    

    #This is experimental part and by default not used by main code
    def ptLikelihood(self, pTTauTau: np.ndarray, type: np.ndarray):

        mask1 = (np.abs(pTTauTau)<0.5)

        if type == 0:
            pT1 = self.leg1P4[:, 0][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 0][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        elif type == 1:
            pT1 = self.leg1P4[:, 1][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 1][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        elif type == 2:
            pT1 = self.leg1P4[:, 2][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
            pT2 = self.leg2P4[:, 2][:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))

        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)[:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))
        x2Min = np.minimum(1.0, self.mVisOverTauSquare2)[:, np.newaxis] * np.ones((1, pTTauTau.shape[1]))

        x1Max = np.ones(pTTauTau.shape)
        x2Max = np.ones(pTTauTau.shape)

        a_x2 = x1Min *pT2/(x1Min*pTTauTau - pT1)
        b_x2 = x1Max*pT2/(x1Max*pTTauTau - pT1)

        x1_singularity = pT1/pTTauTau
        x2_vs_x1_singularity = (x1_singularity>0.0) & (x1_singularity<1.0)

        momentum_sign = (-pT2*pT1<0)

        x2Min = np.where(momentum_sign, np.maximum(x2Min, b_x2), x2Min)
        x2Max = np.where(momentum_sign, np.minimum(x2Max, a_x2), x2Max)
        x2Max = np.where((momentum_sign) & (x2_vs_x1_singularity) & (x2Max<0), 1.0, x2Max)
        x2Min = np.where(~momentum_sign, np.maximum(x2Min, a_x2), x2Min)
        x2Max = np.where(~momentum_sign, np.minimum(x2Max, b_x2), x2Max)
        x2Max = np.where((~momentum_sign) & (x2_vs_x1_singularity) & (x2Max<0), 1.0, x2Max)

        x2Min[x2Min<0] = 0.0
        
        mask2 = (x2Min > x2Max)

        HadDecay1 = np.broadcast_to((self.leg1DecayType != 1)[:, np.newaxis], pTTauTau.shape)
        HadDecay2 = np.broadcast_to((self.leg2DecayType != 1)[:, np.newaxis], pTTauTau.shape)
        
        mNuNuIntegral = np.zeros((pTTauTau.shape))
        x2 = np.minimum(1.0, x2Max)

        term1 = pT2 - pTTauTau*x2
        log_term1 = np.log(np.abs(term1))

        integralMax = pT1*(pTTauTau*x2 + pT2**2/term1 + 2*pT2*log_term1)/pTTauTau**3

        ###MOST CONSUMING PART 1###

        mNuNuIntegral += HadDecay1 * (-pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term1**2 + 6*pT2*log_term1)/(2*pTTauTau**4))
        mNuNuIntegral += HadDecay2 * (-pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term1**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term1 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term1))

        integralMax += mNuNuIntegral

        ###END OF MOST CONSUMING PART 1###

        mNuNuIntegral = np.zeros((pTTauTau.shape))

        x2 = x2Min
        term2 = pT2 - pTTauTau*x2
        log_term2 = np.log(np.abs(term2))

        integralMin = pT1*(pTTauTau*x2+pT2**2/term2+2*pT2*log_term2)/pTTauTau**3

        ###MOST CONSUMING PART 2###
        
        mNuNuIntegral += HadDecay1 * (-pT1**2*(2*pTTauTau*x2+pT2**2*(5*pT2-6*pTTauTau*x2)/term2**2+6*pT2*log_term2)/(2*pTTauTau**4))
        mNuNuIntegral += HadDecay2 * (-pT1/(2*pTTauTau**5)*(2*pT2*pTTauTau*(-3*pT1 + 2*pTTauTau)*x2 + pTTauTau**2*(-pT1 + pTTauTau)*x2**2 + (pT2**4*pT1)/term2**2 + 2*pT2**3*(-4*pT1 + pTTauTau)/term2 + 6*pT2**2*(-2*pT1 + pTTauTau)*log_term2))
        
        integralMin += mNuNuIntegral

        ###END OF MOST CONSUMING PART 2###

        value = integralMax - integralMin

        value[mask1 | mask2] = 0.0

        #value*=1E4

        return np.abs(value)
    
    def metTF(self, metP4: np.ndarray, nuP4: np.ndarray, covMET: np.ndarray):
        aMETx = metP4[..., 0]
        aMETy = metP4[..., 1]

        covDET = np.linalg.det(covMET)
        mask = covDET < 1E-10
        covDET[mask] = 1.0

        constMET = 1/2/np.pi/np.sqrt(covDET)
        residualX = aMETx[:, np.newaxis] - nuP4[:, :, 0]
        residualY = aMETy[:, np.newaxis] - nuP4[:, :, 1]

        #covMET 0 coordinate responds to X and 1 coordinate to Y
        pull2 = residualX*(covMET[:, np.newaxis, 1, 1]*residualX - covMET[:, np.newaxis, 0, 1]*residualY) + residualY*(-covMET[:, np.newaxis, 1, 0]*residualX + covMET[:, np.newaxis, 0, 0]*residualY)
        pull2 /= covDET[:, np.newaxis]
        
        pull2[np.broadcast_to(mask[:, np.newaxis], pull2.shape)] = 0.0
        return constMET[:, np.newaxis]*np.exp(-0.5*pull2)
    

    def value(self, x: np.ndarray):
        
        x1Min = np.minimum(1.0, self.mVisOverTauSquare1)
        x2Min = np.minimum(1.0, self.mVisOverTauSquare2)

        mask = (x[:, 0] < x1Min[:, np.newaxis]) | (x[:, 1] < x2Min[:, np.newaxis])
        
        testP4 = self.leg1P4[:, np.newaxis, :] / x[:, 0][:, np.newaxis] + self.leg2P4[:, np.newaxis, :] / x[:, 1][:, np.newaxis]

        testMET = testP4 - self.leg1P4[:, np.newaxis, :] - self.leg2P4[:, np.newaxis, :]

        value = np.full(testMET.shape[:2], -1.0)
        
        if self.enable_MET:
            value *= self.metTF(self.recoMET, testMET, self.covMET)
        if self.enable_mass:
            value *= self.massLikelihood(InvariantMass(testP4))
        if self.enable_px:
            value *= self.ptLikelihood(testP4[:, :, 0], 0)
        if self.enable_py:
            value *= self.ptLikelihood(testP4[:, :, 1], 1)

        value[mask] = 0.0

        return value

class FastMTT(Likelihood):
    def __init__(self):
        self.myLikelihood = Likelihood()
        self.BestLikelihood = 0.0
        self.BestX = np.array([0.0, 0.0])
        self.bestP4 = 0.0
        self.tau1P4 = 0.0
        self.tau2P4 = 0.0

        #Likelihood of which event to plot
        #-1 = no plot
        self.GridLikelihood = -1

        return
    
    def run(self, measuredTauLeptons, measuredMETx, measuredMETy, covMET) -> np.ndarray:

        start_real_time = time.time()
        start_cpu_time = time.process_time()

        ##############################################
                            #RUN
        ##############################################

        if np.shape(measuredTauLeptons)[1] != 2:
            print(f"Number of MeasuredTauLepton is {len(measuredTauLeptons)}. A user shouls pass exactly two leptons.\n")
            return
        
        metLenght = np.sqrt(measuredMETx**2 + measuredMETy**2)
        aMET = np.array([measuredMETx, measuredMETy, np.zeros(np.shape(measuredMETx)), metLenght]).T


        aLepton1 = measuredTauLeptons[:, 0]
        aLepton2 = measuredTauLeptons[:, 1]
        p4_Lepton1 = self.get_p4(aLepton1)
        p4_Lepton2 = self.get_p4(aLepton2)


        self.myLikelihood.mvisleg1 = aLepton1[:, 4]
        self.myLikelihood.mvisleg2 = aLepton2[:, 4]

        #setMETinputs
        self.myLikelihood.recoMET = aMET
        self.myLikelihood.covMET = covMET

        self.myLikelihood.setLeptonInputs(p4_Lepton1, p4_Lepton2, aLepton1[:, 0], aLepton2[:, 0], aLepton1[:, 5], aLepton2[:, 5])

        self.scan()
        
        self.tau1P4 = p4_Lepton1*(1/self.BestX[:, np.newaxis, 0])
        self.tau2P4 = p4_Lepton2*(1/self.BestX[:, np.newaxis, 1])
        self.bestP4 = self.tau1P4 + self.tau2P4

        ##############################################

        #Time calculation part:
        end_real_time = time.time()
        end_cpu_time = time.process_time()
        
        real_time_elapsed = end_real_time - start_real_time
        cpu_time_elapsed = end_cpu_time - start_cpu_time

        print(f"Real time elapsed: {real_time_elapsed} seconds")
        print(f"CPU time elapsed: {cpu_time_elapsed} seconds")
    
    #lepton[0]: decay_type:
    #1 - TauToHad
    #2 - TauToElec
    #3 - TauToMu

    #lepton[1]: pt
    #lepton[2]: eta
    #lepton[3]: phi
    #lepton[4]: mass
    #lepton[5]: decay_mode (for hadrons):
    #uzupełnić o listę decay modes

    def get_p4(self, lepton: np.ndarray):
        p = lepton[:, 1] * np.cosh(lepton[:, 2])
        px = lepton[:, 1] * np.cos(lepton[:, 3])
        py = lepton[:, 1] * np.sin(lepton[:, 3])
        pz = lepton[:, 1] * np.sinh(lepton[:, 2])
        energy = np.sqrt(p**2 + lepton[:, 4])
        return np.array([px, py, pz, energy]).T
    
    def scan(self):
        
        nGridPoints = 100
        gridFactor = 1.0/nGridPoints

        X1 = np.arange(1, nGridPoints) * gridFactor
        X2 = np.arange(1, nGridPoints) * gridFactor

        # Cartesian product
        pairs = np.column_stack((np.repeat(X1, len(X2)), np.tile(X2, len(X1))))
        
        lh = self.myLikelihood.value(pairs)

        minimum = np.argmin(lh, axis=1)

        self.BestX = pairs[minimum] 
        self.BestLikelihood = lh[np.arange(lh.shape[0]), minimum]

        if self.GridLikelihood != -1:
            self.plot_likelihood(lh, X1, X2, self.GridLikelihood)

        #Code for minimalizing function with scipy.

        '''initial_guess = np.array([0.5, 0.5])
        result = minimize(self.myLikelihood.value, initial_guess, method='BFGS')
        self.BestX = result.x
        self.BestLikelihood = result.fun'''

        
        #Faster than grid search in pure python
        #Slower than grid search in numpy with vectorization and broadcasting
        #Plan to replace it with jax and/or numba

        return
    
    def plot_likelihood(self, lh, X1, X2, event_number=0):
        
        nGridPoints = np.shape(X1)[0]

        lh_grid = lh[event_number, :].reshape(nGridPoints, nGridPoints)

        plt.figure(figsize=(8, 6))
        plt.imshow(lh_grid, origin='lower', extent=(X1.min(), X1.max(), X2.min(), X2.max()), cmap='viridis')
        plt.colorbar(label='Likelihood')
        plt.xlabel("X1")
        plt.ylabel("X2")
        plt.title("2D Heatmap of Likelihood Function")
        plt.savefig("likelihood_heatmap.png", dpi=300)
        plt.close()