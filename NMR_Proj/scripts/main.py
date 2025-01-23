import os
import requests
import pandas as pd
from tqdm import tqdm
from Bio.PDB import PDBParser, Polypeptide
import logging
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Define directories
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # Move one level up
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
BMRB_DIR = os.path.join(DATA_DIR, "bmrb_files")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
os.makedirs(BMRB_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)


def download_bmrb(entry_id):
    """Download BMRB file."""
    url = f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{entry_id}/validation/AVS_full.txt"
    output_path = os.path.join(BMRB_DIR, f"{entry_id}.str")
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            logging.info(f"Downloaded BMRB file for ID {entry_id}")
        else:
            logging.warning(f"BMRB file not found for ID {entry_id}")
            return None
    except Exception as e:
        logging.error(f"Error downloading BMRB file for ID {entry_id}: {e}")
        return None
    return output_path


def download_pdb(pdb_id):
    """Download PDB file."""
    output_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    if os.path.exists(output_path):
        logging.info(f"PDB file for ID {pdb_id} already exists. Skipping download.")
        return output_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            logging.info(f"Downloaded PDB file for ID {pdb_id}")
        else:
            logging.warning(f"PDB file not found for ID {pdb_id}")
            return None
    except Exception as e:
        logging.error(f"Error downloading PDB file for ID {pdb_id}: {e}")
        return None
    return output_path

import re

def parse_chemical_shifts(file_path):
    """
    Parse the chemical shift data for CA, CB, and C values from the AVS file.
    """
    shifts = []
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        current_residue = None

        for line in lines:
            line = line.strip()

            # Check for residue header (e.g., R1, P2)
            if re.match(r"^[A-Z]\d+", line):
                current_residue = {
                    "Residue_Type": re.match(r"^[A-Z]", line).group(),
                    "Residue_ID": int(re.search(r"\d+", line).group()),
                    "C_Shift": None,
                    "CA_Shift": None,
                    "CB_Shift": None,
                }

            # Extract chemical shifts from "Ave C Shift Values" lines
            elif line.startswith("Ave C Shift Values>>"):
                if current_residue:
                    c_shift = re.search(r"C :: ([\d.]+)", line)
                    ca_shift = re.search(r"CA :: ([\d.]+)", line)
                    cb_shift = re.search(r"CB :: ([\d.]+)", line)

                    if c_shift:
                        current_residue["C_Shift"] = float(c_shift.group(1))
                    if ca_shift:
                        current_residue["CA_Shift"] = float(ca_shift.group(1))
                    if cb_shift:
                        current_residue["CB_Shift"] = float(cb_shift.group(1))

                    # Add the current residue to the shifts list
                    shifts.append(current_residue)
                    current_residue = None

        logging.info(f"Parsed {len(shifts)} chemical shifts from {file_path}")
    except Exception as e:
        logging.error(f"Error parsing chemical shifts from {file_path}: {e}")
    return shifts



def parse_secondary_structures(pdb_file):
    """Parse secondary structure data from a PDB file."""
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

        logging.info(f"Parsed secondary structures from {pdb_file}")
    except Exception as e:
        logging.error(f"Error parsing PDB file {pdb_file}: {e}")
    return residues


def save_to_csv(data, output_file):
    """Save data to CSV."""
    try:
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False)
        logging.info(f"Data saved to {output_file}")
    except Exception as e:
        logging.error(f"Error saving data to {output_file}: {e}")


def combine_data(bmrb_file, pdb_file):
    """
    Combine chemical shift and secondary structure data.
    """
    # Parse chemical shifts and secondary structure data
    chem_shifts = parse_chemical_shifts(bmrb_file)  # List of dictionaries
    sec_structures = parse_secondary_structures(pdb_file)  # List of dictionaries

    # Create a dictionary for quick lookup of chemical shifts by Residue_ID
    chem_shifts_dict = {entry["Residue_ID"]: entry for entry in chem_shifts}

    # Combine data
    combined_data = []
    for res in sec_structures:
        res_id = res["Residue_ID"]
        res_type = res["Residue_Type"]
        sec_struct = res["Secondary_Structure"]

        # Retrieve matching chemical shifts
        shifts = chem_shifts_dict.get(res_id, {})
        combined_data.append({
            "Residue_ID": res_id,
            "Residue_Type": res_type,
            "C_Shift": shifts.get("C_Shift"),
            "CA_Shift": shifts.get("CA_Shift"),
            "CB_Shift": shifts.get("CB_Shift"),
            "Secondary_Structure": sec_struct,
        })

    return combined_data


def main():
    bmrb_pdb_list = [
        {"BMRB_ID": "46", "PDB_ID": "1boc"}
    ]

    all_data = []
    for entry in tqdm(bmrb_pdb_list, desc="Processing Entries"):
        bmrb_file = download_bmrb(entry["BMRB_ID"])
        pdb_file = download_pdb(entry["PDB_ID"])

        if not bmrb_file or not pdb_file:
            logging.warning(f"Skipping entry {entry} due to missing files.")
            continue

        combined = combine_data(bmrb_file, pdb_file)
        all_data.extend(combined)

    if all_data:
        output_path = os.path.join(OUTPUT_DIR, "residue_data_table.csv")
        save_to_csv(all_data, output_path)
    else:
        logging.info("No valid data to process or save.")


if __name__ == "__main__":
    main()
