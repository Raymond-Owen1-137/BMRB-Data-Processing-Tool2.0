# BMRB-Data-Processing-Tool2.0
# BMRB Data Processing Tool 2.0

A Python-based tool for automating the extraction, processing, and analysis of chemical shift and secondary structure data from the Biological Magnetic Resonance Bank (BMRB). This tool combines chemical shift values (C, CA, CB) and secondary structure annotations (from PDB files) into a unified format for further analysis.

## Features

- **Automated Downloads**:
  - Retrieves chemical shift data from the BMRB database.
  - Downloads PDB files for secondary structure information.
- **Data Parsing and Integration**:
  - Extracts residue-wise chemical shift values (`C`, `CA`, `CB`).
  - Parses secondary structure data from PDB files.
  - Combines data into a unified CSV file for easy analysis.
- **Error Handling**:
  - Gracefully handles missing or corrupt data files.
- **Output**:
  - Generates clean CSV files with combined data.

---

## Installation

### Prerequisites
Ensure you have the following installed:
- Python 3.7 or later
- `pip` (Python package manager)

### Clone the Repository
```bash
git clone https://github.com/Raymond-Owen1-137/BMRB-Data-Processing-Tool2.0.git
cd BMRB-Data-Processing-Tool2.0
