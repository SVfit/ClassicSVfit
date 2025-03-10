import numpy as np
import uproot
import matplotlib.pyplot as plt
import time
import FastMTT
import multiprocessing as mp

global_fMTT = None  

def init_worker():
    #One FastMTT object for each core
    global global_fMTT
    global_fMTT = FastMTT.FastMTT()

def process_batches_for_worker(worker_batches):
    #Each core processes its own batches
    global global_fMTT
    results = []
    for batch_data in worker_batches:
        measuredTau, METx, METy, covMET = batch_data
        global_fMTT.run(measuredTau, METx, METy, covMET)
        results.append((global_fMTT.mass, global_fMTT.pt, global_fMTT.tau1pt, global_fMTT.tau2pt))
    return results

def process_FastMTT(measuredTauLeptons, xMETs, yMETs, covMETs, batch_size=5_000, num_workers=4):
    num_total = len(measuredTauLeptons)
    num_batches = int(np.ceil(num_total / batch_size))
    
    # Split to cores
    worker_data_splits = np.array_split(range(num_total), num_workers)
    worker_batches = []
    
    for worker_indices in worker_data_splits:
        batches = [
            (measuredTauLeptons[i:i + batch_size],
             xMETs[i:i + batch_size],
             yMETs[i:i + batch_size],
             covMETs[i:i + batch_size])
            for i in range(worker_indices[0], worker_indices[-1] + 1, batch_size)
        ]
        worker_batches.append(batches)
    
    start_time = time.time()
    
    # Multiprocessing
    with mp.Pool(processes=num_workers, initializer=init_worker) as pool:
        results = pool.map(process_batches_for_worker, worker_batches)
    
    # Calculating results
    mFast, ptFast, tau1pt, tau2pt = zip(*[item for sublist in results for item in sublist])
    
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