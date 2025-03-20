import numpy as np
import uproot
from FastMTT_utils import *
import FastMTT

file_path = "/eos/home-w/wmatyszk/public/FastMTT_tutorial.root"
tree_name = "tree;1;1"
branches = ['pt_1', 'eta_1', 'phi_1', 'm_1', 'pt_2', 'eta_2', 'phi_2', 'm_2', 'dm_2', 'met', 'metphi', 'metcov00', 'metcov01', 'metcov11']

data = read_root_file(file_path, tree_name, branches, entry_stop = None)

# Access the numpy arrays
shape = data["pt_1"].shape

measuredTauLeptons = np.array([
    [np.full(shape, 3), data["pt_1"], data["eta_1"], data["phi_1"], data["m_1"], np.full(shape, -1)],
    [np.full(shape, 1), data["pt_2"], data["eta_2"], data["phi_2"], data["m_2"], data["dm_2"]]
])

measuredTauLeptons = np.transpose(measuredTauLeptons, (2, 0, 1))

covMET = np.array([[data["metcov00"], data["metcov01"]], [data["metcov01"], data["metcov11"]]])

covMET = np.transpose(covMET, (2, 0, 1))

METx = data["met"] * np.cos(data["metphi"])
METy = data["met"] * np.sin(data["metphi"])

print("Input shapes: ", measuredTauLeptons.shape, covMET.shape, METx.shape, METy.shape)

fMTT = FastMTT.FastMTT()
fMTT.run(measuredTauLeptons, METx, METy, covMET)
mFast = fMTT.mass
ptFast = fMTT.pt
mFast, ptFast, tau1pt, tau2pt = process_FastMTT(measuredTauLeptons, METx, METy, covMET, batch_size = 10, num_workers = 8)
print("Output shape: ", mFast.shape, ptFast.shape)
print("Output means: ", np.mean(mFast), np.mean(ptFast))