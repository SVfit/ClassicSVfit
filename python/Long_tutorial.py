import numpy as np
#import jax.numpy as jnp
import pandas as pd
import matplotlib.pyplot as plt
import FastMTT
import argparse
import os
from scipy.stats import norm
#import dask.array as da

def process_events_csv(measuredTauLeptons, measuredMETx, measuredMETy, covMET, Higgs_mass, Higgs_pt):

    fMTT = FastMTT.FastMTT()

    #You can choose to plot likelihood for one of the events. -1 means no plot.
    fMTT.WhichLikelihoodPlot = -1

    #You can also choose to calculate uncertainties by:
    fMTT.CalculateUncertainties = True

    #You can enable some likelihood components or other constraints in the following way:
    fMTT.myLikelihood.setWindow = [123, 127] #default settings
    fMTT.myLikelihood.enableLikelihoodComponents(window = True)
    fMTT.myLikelihood.setMassConstraint(mean = 125, sigma = 10) #default settings
    fMTT.myLikelihood.enableLikelihoodComponents(window = False, mass_constraint = True)

    #For default set-up:
    fMTT.myLikelihood.enableLikelihoodComponents(MET = True, mass = True, px = False, py = False, mass_constraint = False, window = False)

    print('Input shapes:', measuredTauLeptons.shape, measuredMETx.shape, measuredMETy.shape, covMET.shape)
    fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
    mFast = fMTT.mass
    print(mFast.shape)

    fastMTT_one_sigma = fMTT.one_sigma
    ptFast = fMTT.pt
    print("FastMTT mass mean:", np.mean(mFast))
    print("FastMTT mass 1 sigma:", np.mean(fastMTT_one_sigma))

    over_range = np.sum((mFast > fMTT.myLikelihood.window[1]))
    below_range = np.sum((mFast < fMTT.myLikelihood.window[0]))
    print(f"Number of cases where Higgs is below {fMTT.myLikelihood.window[0]} GeV: {below_range}")
    print(f"Number of cases where Higgs is above {fMTT.myLikelihood.window[1]} GeV: {over_range}")
    
    print('*****\nFastMTT pt comparison\n*****')
    print('FastMTT pt shape:', ptFast.shape)
    print('True pt shape:', Higgs_pt.shape, '\n')
    
    print("FastMTT pt mean:", np.mean(ptFast))
    print('True pt mean', np.mean(Higgs_pt))
    print('Difference mean:', np.mean(np.absolute(ptFast - Higgs_pt)))

    print('FastMTT pt std:', np.std(ptFast))
    print('True pt std:', np.std(Higgs_pt))
    print('Difference std:', np.std(ptFast - Higgs_pt))

    ###Plots and uncertainty test
    Higgs_plot(mFast)
    pT_resolution_plot(ptFast, Higgs_pt)

    uncertainty_test(Higgs_mass, fMTT)

def load_events_csv(csv_data):

    df = pd.read_csv(csv_data, nrows = 1000)

    event_df = df[['H.m', 'H.pt', 'METx', 'METy', 'covXX', 'covXY', 'covYY', 'dm1', 'pt1', 'eta1', 'phi1', 'mass1', 'type1', 'dm2', 'pt2', 'eta2', 'phi2', 'mass2', 'type2']].copy()

    Higgs_mass = event_df.pop('H.m').to_numpy()
    Higgs_pt = event_df.pop('H.pt').to_numpy()
    METx = event_df.pop('METx').to_numpy()
    METy = event_df.pop('METy').to_numpy()
    metcov = event_df[['covXX', 'covXY', 'covXY', 'covYY']].to_numpy()
    event_df.drop(columns=['covXX', 'covXY', 'covYY'], inplace=True)
    metcov = np.reshape(metcov, (len(metcov), 2, 2))

    print('pandas dataframe:\n', event_df)

    events = event_df.to_numpy()
    events = np.reshape(events, (len(events), 2, 6))

    return {"measuredTauLeptons": events, "measuredMETx": METx, "measuredMETy": METy, "covMET": metcov, "Higgs_mass": Higgs_mass, "Higgs_pt": Higgs_pt}

def Higgs_plot(mFast, output_path = 'images/fastMTT/fastMTT_masses.png'):

    ### ADDITIONAL TESTS AND UI ###
    
    ### PLOTTING ###

    plt.tick_params(axis='both', labelsize=14)

    bin_width = 10
    bins = np.arange(30, 300 + bin_width, bin_width)

    # Calculate the mean and standard deviation
    mean_mass = np.mean(mFast)
    std_mass = np.std(mFast)

    plt.figure(figsize=(8, 6))
    plt.hist(mFast, bins=bins, color='blue', alpha=0.7, edgecolor='black')

    # Add vertical line for the mean
    plt.axvline(mean_mass, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_mass:.2f} GeV')

    # Add vertical lines for one standard deviation
    plt.axvline(mean_mass - std_mass, color='orange', linestyle='--', linewidth=2, label=f'1σ: {std_mass:.2f} GeV')
    plt.axvline(mean_mass + std_mass, color='orange', linestyle='--', linewidth=2)

    # Add labels and title
    plt.xlabel('Mass (GeV)', fontsize=14)
    plt.ylabel('Number of Events', fontsize=14)
    plt.title(f'Two tau leptons invariant mass', fontsize=16)

    # Add grid and legend
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)

    # Save and close the plot
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()

def pT_resolution_plot(ptFast, Higgs_pt, output_path = 'images/fastMTT/fastMTT_pt_differences.png'):

    pt_resolution = (ptFast - Higgs_pt)/Higgs_pt

    ptbin_width = 0.1
    ptbins = np.arange(-3, 3 + ptbin_width, ptbin_width)

    # Drawing a histogram
    plt.figure(figsize=(8, 6))
    plt.hist(pt_resolution, bins = ptbins, color='blue', alpha=0.7, edgecolor='black', density=True)
    plt.title(r'Histogram of pT resolution $(pT_\text{reco} - pT_\text{true})$ \ $pT_\text{true}$ ', fontsize=14)
    plt.xlabel('pT resolution', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)

    # Save to file
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Histogram was saved to {output_path}")

def uncertainty_test(Higgs_mass, fMTT):

    deviations = np.absolute((fMTT.mass - Higgs_mass)) /(fMTT.one_sigma+0.001)
    
    xmin, xmax = -10, 10

    outliers = deviations[(deviations < xmin) | (deviations > xmax)]
    deviations = deviations[(deviations >= xmin) & (deviations <= xmax)]
    
    num_outliers = len(outliers)
    print(f"Number of outliers (>10 sigma): {num_outliers}")

    #chi_square = np.sum(deviations**2) / np.size(deviations)

    within_1sigma = np.sum(deviations <= 1) / np.size(deviations)
    within_2sigma = np.sum(deviations <= 2) / np.size(deviations)
    within_3sigma = np.sum(deviations <= 3) / np.size(deviations)

    print(f"Masses in 1σ: {within_1sigma*100}%")
    print(f"Masses in 2σ: {within_2sigma*100}%")
    print(f"Masses in 3σ: {within_3sigma*100}%")
    #print(f"Chi^2 test: {chi_square}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processes data from a CSV file and prints the results.")
    parser.add_argument("file_path", type=str, help="Path to the CSV file.")
    args = parser.parse_args()

    csv_data = load_events_csv(args.file_path)
    process_events_csv(**csv_data)