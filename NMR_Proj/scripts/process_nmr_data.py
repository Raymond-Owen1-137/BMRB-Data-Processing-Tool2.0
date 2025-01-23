import os
import requests
import pandas as pd
from tqdm import tqdm
from Bio.PDB import PDBParser, Polypeptide
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # Move one level up
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
BMRB_DIR = os.path.join(DATA_DIR, "bmrb_files")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

logging.info("Attempting to download BMRB file...")
logging.error("An error occurred...")

os.makedirs(BMRB_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

def download_bmrb(entry_id):
    print(f"Attempting to download BMRB file for ID {entry_id}...")
    url = f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{entry_id}/validation/AVS_full.txt"
    output_path = os.path.join(BMRB_DIR, f"{entry_id}.str")
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            print(f"Downloaded BMRB file for ID {entry_id}")
        else:
            print(f"BMRB file not found for ID {entry_id}")
            return None
    except Exception as e:
        print(f"Error downloading BMRB file for ID {entry_id}: {e}")
        return None
    return output_path

def download_pdb(pdb_id):
    output_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    if os.path.exists(output_path):
        print(f"PDB file for ID {pdb_id} already exists. Skipping download.")
        return output_path

    print(f"Attempting to download PDB file for ID {pdb_id}...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            print(f"Downloaded PDB file for ID {pdb_id}")
        else:
            print(f"PDB file not found for ID {pdb_id}")
            return None
    except Exception as e:
        print(f"Error downloading PDB file for ID {pdb_id}: {e}")
        return None
    return output_path

import re

def parse_chemical_shifts(bmrb_file):
    shifts = {}
    try:
        with open(bmrb_file, "r") as f:
            lines = f.readlines()
        if not lines:
            print(f"BMRB file {bmrb_file} is empty. Skipping...")
            return shifts

        in_shift_block = False
        for line in lines:
            if "_Atom_chem_shift." in line:
                in_shift_block = True
                continue

            if in_shift_block:
                if line.strip() == "" or line.strip().startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        residue_number = int(parts[1])
                        atom_name = parts[2]
                        shift_value = float(parts[3])
                        if residue_number not in shifts:
                            shifts[residue_number] = {}
                        shifts[residue_number][atom_name] = shift_value
                    except ValueError:
                        continue
        print(f"Parsed chemical shifts from {bmrb_file}: {shifts}")
    except Exception as e:
        print(f"Error parsing BMRB file {bmrb_file}: {e}")
    return shifts

def parse_secondary_structures(pdb_file):
    residues = []
    try:
        secondary_structures = {}
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("HELIX"):
                    start = int(line[21:25].strip())
                    end = int(line[33:37].strip())
                    for res_id in range(start, end + 1):
                        secondary_structures[res_id] = "H"
                elif line.startswith("SHEET"):
                    start = int(line[22:26].strip())
                    end = int(line[33:37].strip())
                    for res_id in range(start, end + 1):
                        secondary_structures[res_id] = "E"

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]

        for chain in model:
            for residue in chain.get_residues():
                if Polypeptide.is_aa(residue, standard=True):
                    res_id = residue.id[1]
                    res_name = residue.resname.strip()
                    sec_struct = secondary_structures.get(res_id, "C")
                    residues.append({
                        "Residue_ID": res_id,
                        "Residue_Type": res_name,
                        "Secondary_Structure": sec_struct
                    })

        print(f"Parsed secondary structures from {pdb_file}: {residues}")
    except Exception as e:
        print(f"Error parsing PDB file {pdb_file}: {e}")
    return residues

def save_to_csv(data, output_file):
    print(f"Saving data to {output_file}...")
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

def combine_data(bmrb_file, pdb_file):
    chem_shifts = parse_chemical_shifts(bmrb_file)
    sec_structures = parse_secondary_structures(pdb_file)

    if not chem_shifts:
        print(f"No chemical shifts found in {bmrb_file}.")
    if not sec_structures:
        print(f"No secondary structure data found in {pdb_file}.")

    combined_data = []
    for res in sec_structures:
        res_id = res["Residue_ID"]
        res_type = res["Residue_Type"]
        sec_struct = res["Secondary_Structure"]

        ca_shift = chem_shifts.get(res_id, {}).get("CA", None)
        cb_shift = chem_shifts.get(res_id, {}).get("CB", None)
        c_shift = chem_shifts.get(res_id, {}).get("C", None)

        combined_data.append({
            "Residue_ID": res_id,
            "Residue_Type": res_type,
            "CA_Shift": ca_shift,
            "CB_Shift": cb_shift,
            "C_Shift": c_shift,
            "Secondary_Structure": sec_struct
        })
    if not combined_data:
        print(f"No valid data for entry {bmrb_file} and {pdb_file}. Skipping...")
    else:
        print(f"Combined data for {bmrb_file} and {pdb_file}: {combined_data}")
    return combined_data

def main():
    bmrb_pdb_list = [
        {"BMRB_ID": "46", "PDB_ID": "1boc"}
    ]

    all_data = []

    for entry in tqdm(bmrb_pdb_list, desc="Processing Entries"):
        print(f"Processing entry: {entry}")
        bmrb_file = os.path.join(BMRB_DIR, f"{entry['BMRB_ID']}.str")
        pdb_file = download_pdb(entry["PDB_ID"])

        if not os.path.exists(bmrb_file):
            print(f"Skipping entry {entry} as the BMRB file is missing.")
            continue

        if not pdb_file:
            print(f"Skipping entry {entry} as the PDB file is missing.")
            continue

        combined = combine_data(bmrb_file, pdb_file)
        all_data.extend(combined)

    if all_data:
        output_path = os.path.join(OUTPUT_DIR, "residue_data_table.csv")
        save_to_csv(all_data, output_path)
        print(f"Data successfully saved to {output_path}")
    else:
        print("No valid data to process or save.")

if __name__ == "__main__":
    main()