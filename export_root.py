import uproot
import numpy as np
import json
import argparse

def root_to_json(root_file_path, output_json_path):
    # Opening the ROOT file
    file = uproot.open(root_file_path)
    data_as_dicts = []

    for key in file.keys():
        obj = file[key]

        # Checking for TTree objects
        if hasattr(obj, "arrays"):
            arrays = obj.arrays(library="np")

            num_rows = len(next(iter(arrays.values())))
            tree_data_as_dicts = [
                {key: arrays[key][i].item() if isinstance(arrays[key][i], np.generic) else arrays[key][i] for key in arrays.keys()}
                for i in range(num_rows)
            ]
            data_as_dicts.extend(tree_data_as_dicts)
            print(f"Added from TTree: {key}")

        else:
            print(f"Ignoring {key} (not TTree)")

    # Saving to JSON
    with open(output_json_path, "w") as f:
        json.dump(data_as_dicts, f, indent=4)

    print(f"TTree saved to {output_json_path}")

# If the script is run directly in the terminal, we could use:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts TTree from ROOT file to JSON")
    parser.add_argument("root_file_path", type=str, help="Path to ROOT file")
    parser.add_argument("output_json_path", type=str, help="Path to JSON file")
    args = parser.parse_args()

    root_to_json(args.root_file_path, args.output_json_path)