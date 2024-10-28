import json

def load_json_data(file_path):
    """
    Wczytuje dane z pliku JSON do Pythona.
    
    Parameters:
        file_path (str): Ścieżka do pliku JSON.
    
    Returns:
        dict: Słownik z danymi JSON.
    """
    with open(file_path, "r") as file:
        data = json.load(file)
    return data

# Przykład użycia
json_data = load_json_data("output.json")
print(json_data)