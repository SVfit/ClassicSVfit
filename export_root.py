import uproot
import numpy as np
import json
import argparse

def root_to_json(root_file_path, output_json_path):
    """
    Konwertuje dane z pliku ROOT do formatu JSON, obsługując obiekty typu TTree oraz TH1D.

    Parameters:
        root_file_path (str): Ścieżka do pliku ROOT.
        output_json_path (str): Ścieżka do wynikowego pliku JSON.
    """
    # Otwieramy plik .root
    file = uproot.open(root_file_path)

    # Inicjalizujemy pustą listę do przechowywania danych
    data_as_dicts = []

    # Iterujemy przez obiekty w pliku
    for key in file.keys():
        obj = file[key]

        # Sprawdzamy, czy obiekt jest typu TTree
        if obj.classname.startswith("TTree"):
            arrays = obj.arrays(library="np")
            num_rows = len(next(iter(arrays.values())))
            data_as_dicts.extend([
                {key: arrays[key][i].item() if isinstance(arrays[key][i], np.generic) else arrays[key][i] for key in arrays.keys()}
                for i in range(num_rows)
            ])

        # Sprawdzamy, czy obiekt jest typu TH1 (histogram)
        elif obj.classname.startswith("TH1"):
            data_as_dicts.append({
                "name": key,
                "bins": obj.values().tolist(),
                "edges": obj.axis().edges().tolist()
            })

    # Zapis do pliku JSON
    with open(output_json_path, 'w') as f:
        json.dump(data_as_dicts, f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Konwertuje plik ROOT do formatu JSON.")
    parser.add_argument("root_file_path", type=str, help="Ścieżka do pliku ROOT.")
    parser.add_argument("output_json_path", type=str, help="Ścieżka do wynikowego pliku JSON.")
    args = parser.parse_args()
    root_to_json(args.root_file_path, args.output_json_path)