import json
import numpy as np
import pandas as pd
import FastMTT
import argparse

def load_json_data(file_path):

    with open(file_path, "r") as file:
        data = json.load(file)
    return data

def process_event_json(json_data):

    for event in json_data:
        
        fMTT = FastMTT.FastMTT()

        measuredMETx = event['metx']
        measuredMETy = event['mety']
        covMET = np.array([[event['metcov00'], event['metcov01']], [event['metcov01'], event['metcov11']]])

        #1 - TauToHad
        #2 - TauToElec
        #3 - TauToMu

        leg1 = np.array([2, event['pt_1'], event['eta_1'], event['phi_1'], event['m_1'], -1])
        leg2 = np.array([1, event['pt_2'], event['eta_2'], event['phi_2'], event['m_2'], event['dm_2']])
        measuredTauLeptons = np.array([leg1, leg2])

        N=100

        measuredTauLeptons = np.tile(measuredTauLeptons, (N, 1, 1))
        measuredMETx = np.tile(measuredMETx, (N, 1))[..., 0]
        measuredMETy = np.tile(measuredMETy, (N, 1))[..., 0]
        covMET = np.tile(covMET, (N, 1, 1))

        fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
        p4Fast = (fMTT.tau1P4 + fMTT.tau2P4)
        mFast = FastMTT.InvariantMass(p4Fast)
        print(f"Fast Mass: {mFast[0]}")

def load_events_csv(csv_data):

    df = pd.read_csv(csv_data)

    event_df = df[['met', 'metphi', 'metcov00', 'metcov01', 'metcov11', 'pt_1', 'eta_1', 'phi_1', 'm_1', 'pt_2', 'eta_2', 'phi_2', 'm_2', 'dm_2']].copy()

    met = event_df.pop('met').to_numpy()
    metphi = event_df.pop('metphi').to_numpy()
    metcov = event_df[['metcov00', 'metcov01', 'metcov01', 'metcov11']].to_numpy()
    event_df.drop(columns=['metcov00', 'metcov01', 'metcov11'], inplace=True)
    metcov = np.reshape(metcov, (len(metcov), 2, 2))

    event_df.insert(0, 'type1', 3)
    event_df.insert(5, 'dm1', -1)
    event_df.insert(6, 'type2', 1)

    events = event_df.to_numpy()
    events = np.reshape(events, (len(events), 2, 6))

    return {"measuredTauLeptons": events, "MET": met, "phiMET": metphi, "covMET": metcov}

def process_events_csv(measuredTauLeptons, MET, phiMET, covMET):

    fMTT = FastMTT.FastMTT()
    
    measuredMETx = MET * np.cos(phiMET)
    measuredMETy = MET * np.sin(phiMET)

    fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
    p4Fast = (fMTT.tau1P4 + fMTT.tau2P4)
    mFast = FastMTT.InvariantMass(p4Fast)

    for i in range(len(mFast)):
        print(f"Fast Mass {i}: {mFast[i]}")

    #Comparison to the C++ results:
    '''df = pd.read_csv('testing_files/results.csv')
    fmtt_df = df['fastMTT_mass'].to_numpy()
    difference = np.absolute(mFast - fmtt_df)
    indices = np.where(difference > 1)[0]
    for index in indices:
        print(f"Event {index}: difference = {difference[index]}")'''
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processes data from either a JSON or CSV file and prints the results.")
    parser.add_argument("file_type", choices=["json", "csv"], help="Type of file to process: 'json' or 'csv'.")
    parser.add_argument("file_path", type=str, help="Path to the file.")
    args = parser.parse_args()

    if args.file_type == "json":
        json_data = load_json_data(args.file_path)
        process_event_json(json_data)
    elif args.file_type == "csv":
        csv_data = load_events_csv(args.file_path)
        process_events_csv(**csv_data)