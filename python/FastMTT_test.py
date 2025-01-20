import numpy as np
#import jax.numpy as jnp
import pandas as pd
import matplotlib.pyplot as plt
import FastMTT
import argparse
import os
from scipy.stats import norm

def load_events_csv(csv_data):

    df = pd.read_csv(csv_data, nrows = 400)

    event_df = df[['H.m', 'METx', 'METy', 'covXX', 'covXY', 'covYY', 'dm1', 'pt1', 'eta1', 'phi1', 'mass1', 'type1', 'dm2', 'pt2', 'eta2', 'phi2', 'mass2', 'type2']].copy()

    Higgs_mass = event_df.pop('H.m').to_numpy()
    METx = event_df.pop('METx').to_numpy()
    METy = event_df.pop('METy').to_numpy()
    metcov = event_df[['covXX', 'covXY', 'covXY', 'covYY']].to_numpy()
    event_df.drop(columns=['covXX', 'covXY', 'covYY'], inplace=True)
    metcov = np.reshape(metcov, (len(metcov), 2, 2))

    print('pandas dataframe:\n', event_df)

    events = event_df.to_numpy()
    events = np.reshape(events, (len(events), 2, 6))

    return {"measuredTauLeptons": events, "measuredMETx": METx, "measuredMETy": METy, "covMET": metcov}

def process_events_csv(measuredTauLeptons, measuredMETx, measuredMETy, covMET):

    fMTT = FastMTT.FastMTT()

    #You can choose to plot likelihood for one of the events. -1 means no plot.
    fMTT.WhichLikelihoodPlot = -1

    #You can also choose to calculate uncertainties by:
    fMTT.CalculateUncertainties = True

    print('Input shapes:', measuredTauLeptons.shape, measuredMETx.shape, measuredMETy.shape, covMET.shape)
    fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
    mFast = fMTT.mass
    fastMTT_one_sigma = fMTT.one_sigma
    print("FastMTT mass mean:", np.mean(mFast))
    print("FastMTT mass 1 sigma:", np.mean(fastMTT_one_sigma))


    ################################################################


    ### ADDITIONAL TESTS AND UI ###



    ### PLOTTING ###

    plt.tick_params(axis='both', labelsize=14)

    bin_width = 10
    bins = np.arange(0, 300 + bin_width, bin_width)

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
    file_path = f"images/fastMTT/fastMTT_histogram.png"
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, dpi=300)
    plt.close()

    ### TEST FOR UNCERTAINTY EVALUATION ###

    uncertainty_test = False
    
    if uncertainty_test:

        deviations = (fMTT.mass - 125) /(fMTT.one_sigma+0.001)
        
        xmin, xmax = -10, 10

        outliers = deviations[(deviations < xmin) | (deviations > xmax)]
        deviations = deviations[(deviations >= xmin) & (deviations <= xmax)]
        
        num_outliers = len(outliers)
        print(f"Number of outliers: {num_outliers}")

        mean, std = norm.fit(deviations)
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mean, std)
        plt.plot(x, p, 'k', linewidth=2, label=f'Normal fit (mean={mean:.2f}, std={std:.2f})')

        plt.hist(deviations, bins=30, density=True, range=(xmin, xmax), alpha=0.6, label="Deviations")
        #x = np.linspace(-100, 100, 1000)
        #plt.plot(x, norm.pdf(x, loc=0, scale=1), 'r-', label="Normal Distribution")
        plt.xlim(xmin, xmax)
        plt.xlabel("Deviation (sigma units)")
        plt.ylabel("Density")
        plt.legend()
        plt.savefig("images/fastMTT/deviations_histogram.png")
        plt.close()

        chi_square = np.sum(deviations**2) / np.size(deviations)

        within_1sigma = np.sum(deviations <= 1) / np.size(deviations)
        within_2sigma = np.sum(deviations <= 2) / np.size(deviations)
        within_3sigma = np.sum(deviations <= 3) / np.size(deviations)

        print(f"Masses in 1σ: {within_1sigma*100}%")
        print(f"Masses in 2σ: {within_2sigma*100}%")
        print(f"Masses in 3σ: {within_3sigma*100}%")
        print(f"Chi^2 test: {chi_square}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processes data from a CSV file and prints the results.")
    parser.add_argument("file_path", type=str, help="Path to the CSV file.")
    args = parser.parse_args()

    csv_data = load_events_csv(args.file_path)
    process_events_csv(**csv_data)