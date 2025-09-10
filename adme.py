import warnings # ignore warnings
warnings.filterwarnings("ignore", category = RuntimeWarning)

from adme_py import ADME # import adme-py package and specifically ADME
import pubchempy as pcp # import PubChem package
import pandas as pd # import Pandas library
import os

file_name = input("Enter the spreadsheet file name: ")
df = pd.read_excel(file_name)

df["Canonical SMILES"] = "" # make empty keys in dataframe dictionary
df["Lipinski's Rule of Five"] = ""
df["BBB Permeable"] = ""
df["Final Decision"] = ""

for index, row in df.iterrows():
    compound_name = row["Compound"] # get the compound name whose adme is to be checked
    compound = pcp.get_compounds(compound_name, 'name') # check on PubChem for data
    if not compound: # if the compound is not found, print an error message
        df.at[index, "Canonical SMILES"] = "NOT FOUND!" # add to dataframe
        continue
    else: # if the compound is found, do the rest
        ligand = compound[0] # get the zeroeth element of the list
        smiles = ligand.connectivity_smiles # get the canonical smiles
        
        summary = ADME(smiles).calculate() # extract and store the ADME results of the particular molecule represented as SMILES in the dictionary summary
        
        lipinski = summary['druglikeness']['lipinski'] # to determine whether the particular ligand passes the lipinski's rules or not
        bbb_permeable = summary['pharmacokinetics']['blood_brain_barrier_permeant'] # to determine whether the particular ligand is blood brain barrier permeant or not
        
        if str(lipinski) == "Pass" and str(bbb_permeable) == "False": # selected or not
            selected = "PERFECT!"
        elif str(lipinski) != "Pass" and str(bbb_permeable) == "False":
            selected = f"Violates {len(lipinski.keys())} Lipinski's Rule(s)!"
        elif str(lipinski) == "Pass" and str(bbb_permeable) == "True":
            selected = "Blood-Brain Barrier Permeable!"
        else:
            selected = "NO"

        df.at[index, "Canonical SMILES"] = smiles # store all the data in the dataframe
        df.at[index, "Lipinski's Rule of Five"] = lipinski
        df.at[index, "BBB Permeable"] = bbb_permeable
        df.at[index, "Final Decision"] = selected
        
output_file = file_name.replace(".xlsx", "") + "_results.xlsx"
df.to_excel(output_file, index = False)

if os.path.exists(output_file) and os.path.getsize(output_file) != 0:
    print(f"{output_file} has been created and saved successfully!")