import numpy as np
import uproot
import matplotlib.pyplot as plt
import time
import FastMTT
import multiprocessing as mp

global_fMTT = None  

def init_worker():
    global global_fMTT
    global_fMTT = FastMTT.FastMTT()

def process_batch(batch_data):
    global global_fMTT  # Each process uses its own class instance

    measuredTau, METx, METy, covMET = batch_data
    global_fMTT.run(measuredTau, METx, METy, covMET)
    return global_fMTT.mass, global_fMTT.pt, global_fMTT.tau1pt, global_fMTT.tau2pt

def process_FastMTT(measuredTauLeptons, xMETs, yMETs, covMETs, batch_size=125, num_workers=4):

    # Preparing batches
    num_batches = int(np.ceil(len(measuredTauLeptons) / batch_size))
    batches = [
        (measuredTauLeptons[i * batch_size:(i + 1) * batch_size],
         xMETs[i * batch_size:(i + 1) * batch_size],
         yMETs[i * batch_size:(i + 1) * batch_size],
         covMETs[i * batch_size:(i + 1) * batch_size])
        for i in range(num_batches)
    ]

    start_time = time.time()

    # Pool multiprocessing with FastMTT instance initialization
    with mp.Pool(processes=num_workers, initializer=init_worker) as pool:
        results = pool.map(process_batch, batches)

    # Ziping results
    mFast, ptFast, tau1pt, tau2pt = zip(*results)
    
    end_time = time.time()
    print(f"Processing FastMTT took {end_time - start_time:.2f} seconds")

    return np.concatenate(mFast, axis=0), np.concatenate(ptFast, axis=0), np.concatenate(tau1pt, axis=0), np.concatenate(tau2pt, axis=0)

def read_root_file(file_path, tree_name, branches, entry_stop=None):
    # Open the ROOT file using uproot
    with uproot.open(file_path) as file:
        # Get the tree from the file
        tree = file[tree_name]
        
        # Read the branches into numpy arrays
        data = tree.arrays(branches, library="np", entry_stop=entry_stop)
    
    return data