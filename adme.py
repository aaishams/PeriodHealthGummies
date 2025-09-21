import warnings # ignore warnings
warnings.filterwarnings("ignore", category = RuntimeWarning)

from adme_py import ADME # import adme-py package and specifically ADME
import pubchempy as pcp # import PubChem package
import pandas as pd # import Pandas library
import os

def final_decision(file, out):
    df = pd.read_excel(file)

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
    df.to_excel(out, index = False)

def detailed_output(file, details):
    df2 = pd.read_excel(file)
    with open(details, "w") as opened:
        opened.write("Detailed Summary of ADME Studies: \n")
        opened.write("\n")

        for index, row in df2.iterrows():
            compound_name = row["Compound"] # get the compound name whose adme is to be checked
            compound = pcp.get_compounds(compound_name, 'name') # check on PubChem for data
            print(index + 1, compound_name)
            opened.writelines(str(index + 1) + ". ")
            opened.writelines(compound_name)
            opened.write("\n")
            if not compound:
                opened.write("NOT FOUND! \n")
                opened.write("\n")
                continue
            else:
                ligand = compound[0]
                smiles = ligand.connectivity_smiles

                summary = ADME(smiles).calculate() # extract and store the ADME results of the particular molecule represented as SMILES in the dictionary summary
                print(summary)
                for key, value in summary.items():
                    opened.write(f"{key}: {value} \n")
                opened.write("\n")

file_name = input("Enter the spreadsheet file name: ")
output_file = file_name.replace(".xlsx", "") + "_results.xlsx"
detailed_file = file_name.replace(".xlsx", "") + "_details.txt"
df = final_decision(file_name, output_file)
df2 = detailed_output(file_name, detailed_file)

if os.path.exists(output_file) and os.path.getsize(output_file) != 0:
    print(f"{output_file} has been created and saved successfully!")
if os.path.exists(detailed_file) and os.path.getsize(detailed_file) != 0:
    print(f"{detailed_file} has been created and saved successfully!")