###PROTOTYPE###
#This is the test code to check if python port works properly
#For now it is only for MeasuredTauLepton class
#And it seems the class works fine

#The next step will be to rewrite it to the test code for full FastMTT

import json
import numpy as np
import vector
from MeasuredTauLepton import MeasuredTauLepton
import argparse

def load_json_data(file_path):

    with open(file_path, "r") as file:
        data = json.load(file)
    return data

def process_events(json_data):
    
    i = 0

    for event in json_data:
        if i>5:
            break

        i += 1

        phiMET = event['metphi']
        measuredMETx = event['met'] * np.cos(phiMET)
        measuredMETy = event['met'] * np.sin(phiMET)
        
        covMET = np.array([
            [event['metcov00'], event['metcov01']],
            [event['metcov01'], event['metcov11']]
        ])
        
        measuredTauLeptons = [
            MeasuredTauLepton(3, event['pt_1'], event['eta_1'], event['phi_1'], event['m_1'], -1),  # The last integer is a decay mode
            MeasuredTauLepton(1, event['pt_2'], event['eta_2'], event['phi_2'], event['m_2'], event['dm_2'])
        ]
        
        print(f"Processed event with run: {event['run']}, lumi: {event['lumi']}, evt: {event['evt']}")
        print(f"Type of first event: {measuredTauLeptons[0].type}")
        print(f"phi_2 of second event: {measuredTauLeptons[1].phi}")
        print(f"measuredMETx: {measuredMETx}, measuredMETy: {measuredMETy}")
        print(f"Determinant of covMET: {np.linalg.det(covMET)}")


# If the script is run directly in the terminal, we could use:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Przetwarza dane z pliku JSON i wypisuje wyniki.")
    parser.add_argument("json_file_path", type=str, help="Ścieżka do pliku JSON.")
    args = parser.parse_args()

    json_data = load_json_data(args.json_file_path)

    process_events(json_data)