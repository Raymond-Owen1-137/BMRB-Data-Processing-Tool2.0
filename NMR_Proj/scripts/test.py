import pandas as pd

# Data as a Python list of dictionaries
data = [
    {"Residue_ID": 0, "Residue_Type": "MET", "C_Shift": None, "CA_Shift": None, "CB_Shift": None, "Secondary_Structure": "C"},
    {"Residue_ID": 1, "Residue_Type": "LYS", "C_Shift": 176.49, "CA_Shift": 56.95, "CB_Shift": 30.66, "Secondary_Structure": "C"},
    # ... (the rest of the data)
]

# Convert to DataFrame and save
df = pd.DataFrame(data)
df.to_csv("residue_data_table.csv", index=False)
print("Data saved to residue_data_table.csv")

filtered_data = df.dropna(subset=["C_Shift", "CA_Shift", "CB_Shift"])
print(filtered_data)
