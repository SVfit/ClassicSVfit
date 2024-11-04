import json
import numpy as np
import vector
import MeasuredTauLepton as mtl
import FastMTT
import argparse

def load_json_data(file_path):

    with open(file_path, "r") as file:
        data = json.load(file)
    return data

def process_event(json_data):
    i = 0
    for event in json_data:
        if i > 5:
            break
        
        fMTT = FastMTT.FastMTT()

        i += 1
        measuredMETx = event['metx']
        measuredMETy = event['mety']
        covMET = np.array([[event['metcov00'], event['metcov01']], [event['metcov01'], event['metcov11']]])

        #1 - TauToHad
        #2 - TauToElec
        #3 - TauToMu

        leg1 = mtl.MeasuredTauLepton(2, event['pt_1'], event['eta_1'], event['phi_1'], event['m_1'], -1) 
        leg2 = mtl.MeasuredTauLepton(1, event['pt_2'], event['eta_2'], event['phi_2'], event['m_2'], event['dm_2'])
        measuredTauLeptons = [leg1, leg2]

        fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
        p4Fast = (fMTT.tau1P4 + fMTT.tau2P4)
        mFast = p4Fast.mass
        print(f"Fast Mass: {mFast}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processes data from a JSON file and prints the results.")
    parser.add_argument("json_file_path", type=str, help="Path to the JSON file.")
    args = parser.parse_args()

    json_data = load_json_data(args.json_file_path)

    process_event(json_data)