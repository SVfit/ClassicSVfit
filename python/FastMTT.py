#import dask.array as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import time
from scipy.constants import physical_constants
import os
import Likelihood
#from scipy.optimize import minimize

###Main reference: https://github.com/SVfit/ClassicSVfit/blob/fastMTT_2024/src/FastMTT.cc ###


ElectronMass = physical_constants['electron mass energy equivalent in MeV'][0]/1000 #MeV -> GeV
MuonMass = physical_constants['muon mass energy equivalent in MeV'][0]/1000 #MeV -> GeV
ChargedPionMass = 139.5/1000 #MeV -> GeV


#Invariant mass calculation
def InvariantMass(p4):
    metric = np.array([-1,-1,-1,1])
    p4_square = p4*(metric*p4)
    m = np.sqrt(np.sum(p4_square, axis=-1))
    return m

def pT(aP4):
    return np.sqrt(aP4[..., 0]**2 + aP4[..., 1]**2)

class FastMTT:
    def __init__(self, calculate_uncertainties = False):
        self.myLikelihood = Likelihood.Likelihood()
        self.BestLikelihood = 0.0
        self.BestX = np.array([0.0, 0.0])
        self.bestP4 = 0.0
        self.tau1P4 = 0.0
        self.tau2P4 = 0.0
        self.mass = 0.0

        #New component to calculate uncertainty
        #It produces long tails, but apart from that calculates uncertainties event by event quite ok ~ after some cuts results are aprox. Gaussian
        #A bit time consuming -- doubles the time of calculation -- so it is disabled by default

        self.CalculateUncertainties = calculate_uncertainties
        self.one_sigma = 0.0

        #Number of event for which likelihood plot will be shown.
        #-1 = no plot
        self.WhichLikelihoodPlot = -1

        return
    
    def run(self, measuredTauLeptons, measuredMETx, measuredMETy, covMET):

        start_real_time = time.time()
        start_cpu_time = time.process_time()

        ##############################################
                            #RUN
        ##############################################

        if np.shape(measuredTauLeptons)[1] != 2:
            print(f"Number of MeasuredTauLepton is {len(measuredTauLeptons)}. A user should pass exactly two leptons.\n")
            return

        metLenght = np.sqrt(measuredMETx**2 + measuredMETy**2)
        aMET = np.array([measuredMETx, measuredMETy, np.zeros(np.shape(measuredMETx)), metLenght]).T

        aLepton1 = measuredTauLeptons[:, 0]
        aLepton2 = measuredTauLeptons[:, 1]

        self.p4_Lepton1 = self.get_p4(aLepton1)
        self.p4_Lepton2 = self.get_p4(aLepton2)

        self.Lepton1 = self.p4_Lepton1
        self.Lepton2 = self.p4_Lepton2

        aLepton1 = self.modify_lepton_mass(aLepton1)
        aLepton2 = self.modify_lepton_mass(aLepton2)

        self.myLikelihood.mvisleg1 = aLepton1[:, 4]
        self.myLikelihood.mvisleg2 = aLepton2[:, 4]

        #setMETinputs
        self.myLikelihood.recoMET = aMET
        self.myLikelihood.covMET = covMET

        self.myLikelihood.setLeptonInputs(self.p4_Lepton1, self.p4_Lepton2, aLepton1[:, 0], aLepton2[:, 0], aLepton1[:, 5], aLepton2[:, 5])

        self.scan()
        
        self.tau1P4 = self.p4_Lepton1*(1/self.BestX[:, np.newaxis, 0])
        self.tau2P4 = self.p4_Lepton2*(1/self.BestX[:, np.newaxis, 1])

        if self.myLikelihood.enable_window:
            mvis = InvariantMass(self.p4_Lepton1 + self.p4_Lepton2)
            mask = mvis > self.myLikelihood.window[1]
            self.tau1P4[mask] = self.p4_Lepton1[mask]
            self.tau2P4[mask] = self.p4_Lepton2[mask]

        self.bestP4 = self.tau1P4 + self.tau2P4
        self.mass = InvariantMass(self.bestP4)
        self.pt = pT(self.bestP4)

        if self.myLikelihood.enable_window:
            self.mass[(self.mass < self.myLikelihood.window[0])] = self.myLikelihood.window[0]
            self.mass[(self.mass > self.myLikelihood.window[1])] = self.myLikelihood.window[1]

        self.tau1pt = np.sqrt(self.tau1P4[..., 0]**2 + self.tau1P4[..., 1]**2)
        self.tau2pt = np.sqrt(self.tau2P4[..., 0]**2 + self.tau2P4[..., 1]**2)

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
    #lepton[5]: hadron decay mode

    def get_p4(self, lepton):
        p = lepton[:, 1] * np.cosh(lepton[:, 2])
        px = lepton[:, 1] * np.cos(lepton[:, 3])
        py = lepton[:, 1] * np.sin(lepton[:, 3])
        pz = lepton[:, 1] * np.sinh(lepton[:, 2])
        energy = np.sqrt(p**2 + lepton[:, 4])
        return np.array([px, py, pz, energy]).T
    
    def modify_lepton_mass(self, aLepton1, electron_mass=ElectronMass, muon_mass=MuonMass, pion_mass=ChargedPionMass):

        # 1) Set electron inv mass to electron
        aLepton1[aLepton1[:, 0] == 2, 4] = electron_mass

        # 2) Set muon inv mass to muon
        aLepton1[aLepton1[:, 0] == 3, 4] = muon_mass

        # 3) If decay is hadronic, we set minimal and maximal mass, depending on decay mode
        mask_type1 = aLepton1[:, 0] == 1

        # a) If hadronic decay mode == -1, then min mass = pion_mass and max mass = 1.5
        mask_a = (mask_type1 & (aLepton1[:, 5] == -1) & (aLepton1[:, 4] < pion_mass))
        aLepton1[mask_a, 4] = pion_mass

        mask_b = (mask_type1 & (aLepton1[:, 5] == -1) & (aLepton1[:, 4] > 1.5))
        aLepton1[mask_b, 4] = 1.5

        # c) If hadronic decay mode == 0, then mass = pion mass
        mask_c = (mask_type1 & (aLepton1[:, 5] == 0))
        aLepton1[mask_c, 4] = pion_mass

        # d) If hadronic decay mode has another value, we set min mass to 0.3 and max mass to 1.5
        mask_d = (mask_type1 & (aLepton1[:, 5] != -1) & (aLepton1[:, 5] != 0) & (aLepton1[:, 4] < 0.3))
        aLepton1[mask_d, 4] = 0.3

        mask_e = (mask_type1 & (aLepton1[:, 5] != -1) & (aLepton1[:, 5] != 0) & (aLepton1[:, 4] > 1.5))
        aLepton1[mask_e, 4] = 1.5

        return aLepton1
    
    def scan(self):
        
        nGridPoints = 100
        gridFactor = 1.0/nGridPoints

        X1 = np.arange(1, nGridPoints+1) * gridFactor
        X2 = np.arange(1, nGridPoints+1) * gridFactor

        # Cartesian product
        self.pairs = np.column_stack((np.repeat(X1, len(X2)),
                                 np.tile(X2, len(X1))))
        
        self.lh = self.myLikelihood.value(self.pairs)

        minimum = np.argmin(self.lh, axis=1)

        self.BestX = self.pairs[minimum]
        self.BestLikelihood = self.lh[np.arange(self.lh.shape[0]), minimum]

        ### USER INTERFACE AND ADDITIONAL COMPONENTS ###

        chi_square = 2.3

        ###
        # 1 sigma = 2.3
        # 2 sigma = 6.0
        # 3 sigma = 9.2
        # 5 sigma = 28.7
        ###

        #Plotting likelihoods
        if self.WhichLikelihoodPlot != -1:
            threshold=self.BestLikelihood[self.WhichLikelihoodPlot]/np.exp(chi_square/2)
            self.plot_likelihood(X1, X2, event_number = self.WhichLikelihoodPlot, threshold=threshold)

        if self.CalculateUncertainties == True:
            self.contour_uncertainties(X1, X2, chi_square)

        #Code for minimalizing function with scipy:

        '''initial_guess = np.array([0.5, 0.5])
        result = minimize(self.myLikelihood.value, initial_guess, method='BFGS')
        self.BestX = result.x
        self.BestLikelihood = result.fun'''
        
        #Faster than grid search in pure python
        #Slower than grid search in numpy with vectorization and broadcasting
        #Potentially one can replace it with jax and/or numba

        #self.check_likelihood()

        return
    
    def plot_likelihood_subplot(self, ax, X1, X2, lh_grid, maximum_likelihood, threshold=None):

        ###Function used by plot_likelihood###
        #(written in this way to make easier to implement plotting multiple likelihoods)

        # Main heatmap
        img = ax.imshow(-lh_grid, origin='lower', extent=(X1[0], X1[-1], X2[0], X2[-1]), 
                        cmap='viridis_r', interpolation='nearest')

        ax.set_aspect('auto')

        # Excluded regions
        unphysical_region_mask = (lh_grid == 0.000001)
        ax.contourf(X1, X2, unphysical_region_mask, levels=[0.5, 1], colors='red', alpha=0.6)
        contour_proxy = plt.Line2D([0], [0], linestyle="none", marker="s", markersize=10, markerfacecolor="red", alpha=0.6)
        ax.legend([contour_proxy], ["Kinematic limit"], loc="upper left")

        # Max. likelihood point
        ax.scatter(maximum_likelihood[1], maximum_likelihood[0], color='red', marker='x', s=100, label="Maximum likelihood", linewidth=2.5)

        # Contour for uncertainty
        if threshold is not None:
            ax.contour(X1, X2, lh_grid, levels=[threshold], colors='blue', linewidths=2, linestyles='--')

        ax.set_xlabel(r"$X_1$")
        ax.set_ylabel(r"$X_2$")
        ax.set_title("2D Heatmap of Likelihood Function")

        return img


    def plot_likelihood(self, X1, X2, event_number=0, threshold=None):
        print("Threshold: ", threshold)
        nGridPoints = np.shape(X1)[0]
        lh_grid = self.lh[event_number, :].reshape(nGridPoints, nGridPoints)
        maximum_likelihood = self.BestX[event_number]
        
        fig, ax = plt.subplots(figsize=(6, 4.5))
        img = self.plot_likelihood_subplot(ax, X1, X2, lh_grid, maximum_likelihood, threshold)

        fig.colorbar(img, ax=ax, label='Likelihood')

        params = {'legend.fontsize': 'xx-large',
                'figure.figsize': (10, 7),
                'axes.labelsize': 'xx-large',
                'axes.titlesize': 'xx-large',
                'xtick.labelsize': 'xx-large',
                'ytick.labelsize': 'xx-large'}
        plt.rcParams.update(params)

        file_path = f"images/fastMTT/likelihood_{self.WhichLikelihoodPlot}_event.png"
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        plt.savefig(file_path, format='png')
        plt.close()

    def evaluate_mass(self, x):
        tau1P4 = self.p4_Lepton1[:, np.newaxis, np.newaxis, :]*(1/x[np.newaxis, :, :, np.newaxis, 0])
        tau2P4 = self.p4_Lepton2[:, np.newaxis, np.newaxis, :]*(1/x[np.newaxis, :, :, np.newaxis, 1])
        bestP4 = tau1P4 + tau2P4
        mass = InvariantMass(bestP4)
        return mass
    
    def contour_uncertainties(self, X1, X2, chi_square = 2.3):
        threshold = self.BestLikelihood/np.exp(chi_square/2)

        nGridPoints = np.shape(X1)[0]
        nEvents = np.shape(self.lh)[0]

        lh_grid = self.lh.reshape(nEvents, nGridPoints, nGridPoints)
        pairs = self.pairs.reshape(nGridPoints, nGridPoints, 2)
        mask = (lh_grid < threshold[:, np.newaxis, np.newaxis])

        up = np.roll(mask, shift=-1, axis=1)
        down = np.roll(mask, shift=1, axis=1)
        left = np.roll(mask, shift=-1, axis=2)
        right = np.roll(mask, shift=1, axis=2)
        boundary_mask = mask & ~(up & down & left & right)
        
        masses = self.evaluate_mass(pairs)
        masses = np.where(boundary_mask, masses, np.nan)
        self.max_masses = np.nanmax(masses, axis=(1, 2))
        self.min_masses = np.nanmin(masses, axis=(1, 2))

        self.one_sigma = (self.max_masses - self.min_masses)/2