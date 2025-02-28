import numpy as np
import uproot
import matplotlib.pyplot as plt
import time
import FastMTT

def process_FastMTT(measuredTauLeptons, xMETs, yMETs, covMETs, batch_size = 5_000, log_interval = 1):
    
    # Splitting into batches
    measuredTauLeptons_batches = np.array_split(measuredTauLeptons, np.ceil(len(measuredTauLeptons) / batch_size))
    METx_batches = np.array_split(xMETs, np.ceil(len(xMETs) / batch_size))
    METy_batches = np.array_split(yMETs, np.ceil(len(yMETs) / batch_size))
    covMET_batches = np.array_split(covMETs, np.ceil(len(covMETs) / batch_size))

    fMTT = FastMTT.FastMTT()

    mFast = []
    ptFast = []

    start_time = time.time()

    for i, (measuredTau, METx, METy, covMET) in enumerate(zip(measuredTauLeptons_batches, METx_batches, METy_batches, covMET_batches)):
        #Each batch processing
        fMTT.run(measuredTau, METx, METy, covMET)
        mFast.append(fMTT.mass)
        ptFast.append(fMTT.pt)

        if i % log_interval == 0 or i == len(measuredTauLeptons_batches) - 1:
            print(f"Batch {i+1}/{len(measuredTauLeptons_batches)} processed")

    #Time measurement
    end_time = time.time()
    print(f"Processing FastMTT took {end_time - start_time:.2f} seconds")

    # Merging:)
    return np.concatenate(mFast, axis=0), np.concatenate(ptFast, axis=0)
    

def read_root_file(file_path, tree_name, branches, entry_stop=None):
    # Open the ROOT file using uproot
    with uproot.open(file_path) as file:
        # Get the tree from the file
        tree = file[tree_name]
        
        # Read the branches into numpy arrays
        data = tree.arrays(branches, library="np", entry_stop=entry_stop)
    
    return data