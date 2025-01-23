import os
import requests
import pandas as pd
from tqdm import tqdm

# Define the base directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # Move one level up
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
BMRB_DIR = os.path.join(DATA_DIR, "bmrb_files")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")

os.makedirs(BMRB_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

def download_bmrb(entry_id):
    """Download BMRB file for the given entry ID."""
    print(f"Attempting to download BMRB file for ID {entry_id}...")
    url = f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{entry_id}/validation/AVS_full.txt"
    output_path = os.path.join(BMRB_DIR, f"{entry_id}_AVS_full.txt")
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            print(f"Downloaded BMRB file for ID {entry_id}")
        else:
            print(f"BMRB file not found for ID {entry_id} (HTTP {response.status_code})")
            return None
    except requests.RequestException as e:
        print(f"Error downloading BMRB file for ID {entry_id}: {e}")
        return None
    return output_path

def parse_avs_file(file_path):
    """Parse AVS data from a text file."""
    results = []
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        for line in lines:
            # Process lines that are not empty and start with an alphabetic character
            if line.strip() and line[0].isalpha():
                results.append(line.strip())

        print(f"Parsed AVS data from {file_path}")
    except Exception as e:
        print(f"Error parsing AVS file {file_path}: {e}")
    return results

def save_avs_to_file(data, output_file):
    """Save parsed AVS data to a file."""
    try:
        with open(output_file, "w") as f:
            f.write("\n".join(data))
        print(f"AVS data saved to {output_file}")
    except Exception as e:
        print(f"Error saving AVS data to file {output_file}: {e}")

def main():
    """Main function to download and process AVS data."""
    bmrb_id = "46"
    avs_file = download_bmrb(bmrb_id)

    if avs_file:
        parsed_data = parse_avs_file(avs_file)
        output_file = os.path.join(OUTPUT_DIR, "parsed_avs_data.txt")
        save_avs_to_file(parsed_data, output_file)
    else:
        print("No AVS file to process.")

if __name__ == "__main__":
    main()
