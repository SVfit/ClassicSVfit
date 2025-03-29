from pybind_wrapper import *
import ROOT as root


decayType = 2
leps1_pt = 33.7393
leps1_eta = 0.9409
leps1_phi = -0.541458
leps1_mass = 0.51100e-3
leps1_dm = 0

leps2_pt = 25.7322
leps2_eta = 0.618228
leps2_phi = 2.79362
leps2_mass = 0.13957
leps2_dm = 0


measuredMETx =  11.7491
measuredMETy = -51.9172

# define MET covariance
covMET = root.TMatrixD(2, 2)
covMET[0][0] =  787.352
covMET[1][0] = -178.63
covMET[0][1] = -178.63
covMET[1][1] =  179.545


leps1_MeasuredTauLepton = MeasuredTauLepton(decayType, leps1_pt, leps1_eta, leps1_phi ,leps1_mass, leps1_dm) # tau -> e (Pt, eta, phi, mass)
leps2_MeasuredTauLepton = MeasuredTauLepton(1, leps2_pt, leps2_eta, leps2_phi ,leps2_mass, leps2_dm) # tau hadronic decay 0 (Pt, eta, phi, mass)



# FastMTT Algorithm
FastMTTAlgo = FastMTT()


leps_MeasuredTauLepton = [leps1_MeasuredTauLepton, leps2_MeasuredTauLepton]

# Your code inside the for loop
FastMTTAlgo.run(leps_MeasuredTauLepton, measuredMETx,  measuredMETy, covMET)

# Get the tau1 and tau2 P4 from FastMTT as LorentzVector
tau1P4mtt = FastMTTAlgo.getTau1P4()
tau2P4mtt = FastMTTAlgo.getTau2P4()

print(f"tau1P4mtt : {tau1P4mtt}")
print(f"tau1P4mtt : {tau1P4mtt}")


# Get the ditau mass
ditauMass = (tau1P4mtt + tau2P4mtt).M()

print(f"ditauMass : {ditauMass}")
  
