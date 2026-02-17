"""
Script for searching PubChem molecules against SMARTS patterns with OpenEye tools.
This script processes molecular structures from PubChem data files, filters them based on
allowed elements and heavy atom counts, and matches them against SMARTS patterns derived
from nitrogen, sulfur, and phosphorous datasets. For each SMARTS pattern, it collects up to
n_mols molecules that match the pattern and have low structural similarity (Tanimoto < 0.4)
to previously selected molecules.
Key Functions:

    get_tree_fp: Generates a Tree fingerprint for a molecule using OpenEye.
Workflow:
    1. Load SMARTS patterns from CSV files (nitrogen, sulfur, phosphorous).
    2. Read molecules from a PubChem SDF file stream.
    3. Filter molecules by element composition, heavy atom count, and custom filter rules.
    4. For each molecule, attempt to match it against SMARTS patterns.
    5. For matched patterns, add the molecule only if it exhibits low similarity to existing molecules.
    6. Export results to a JSON file mapping SMARTS patterns to selected molecule SMILES strings.
Note on "iterate over a shallow copy":
    The line `for pattern in all_smarts[:]` uses a slice `[:]` to create a shallow copy of the
    `all_smarts` list. This allows safe removal of patterns from the original list during iteration
    (via `all_smarts.remove(pattern)`) without causing iteration errors. Iterating directly over
    `all_smarts` while modifying it would result in skipped elements or unexpected behavior.
"""
import pandas as pd
import json
from openeye import oechem
from collections import defaultdict
from openeye import oegraphsim
from openeye import oemolprop


def get_tree_fp(oemol):
    fp = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeFP(fp, oemol, oegraphsim.OEFPType_Tree)
    return fp


n_df = pd.read_csv('nitrogen_summary_updated.csv', delimiter=',')
n_df = n_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

s_df = pd.read_csv('sulfur_summary_updated.csv', delimiter=',')
s_df = s_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

p_df = pd.read_csv('phosphorous_summary_updated.csv', delimiter=',')
p_df = p_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

n_smarts = n_df['SMARTS'].to_list()
s_smarts = s_df['SMARTS'].to_list()
p_smarts = p_df['SMARTS'].to_list()

all_smarts = p_smarts + s_smarts + n_smarts
n_mols = 100  # Maximum number of molecules to collect per SMARTS pattern, per charge category
fname = 'file_name'
ifs = oechem.oemolistream(f'/dfs6/pub/pbehara/check_chembl_parameter_coverage/pubchem/Compounds/{fname}')

# 1. Create an input stream from your custom filter file
filter_file_path = "FILTER.OE"
filter_istream = oechem.oeifstream(filter_file_path)

# Check if the file was successfully opened
if not filter_istream.IsValid():
    print(f"Error: Could not open custom filter file {filter_file_path}")
    exit()

# 2. Initialize the OEFilter object with the custom filter stream
filter_obj = oemolprop.OEFilter(filter_istream)
print(filter_obj.GetTypeCheck())
filter_obj.SetTypeCheck(True)

# Pre-compile OESubSearch objects for all SMARTS patterns
smarts_search_dict = {}
for pattern in all_smarts:
    ss = oechem.OESubSearch(pattern)
    smarts_search_dict[pattern] = ss

mol = oechem.OEGraphMol()
smarts_dict = defaultdict(dict)
mols_added = []
while oechem.OEReadMolecule(ifs, mol):
    added_to_pattern = False
    oechem.OEDeleteEverythingExceptTheFirstLargestComponent(mol)
    if not filter_obj(mol): 
        continue
    mol_name = mol.GetTitle()
    smiles = oechem.OEMolToSmiles(mol)
    if smiles in mols_added:
        continue
    for pattern in all_smarts[:]:  # Iterate over a shallow copy
        # skip already added molecules
        if added_to_pattern:
            break
        if len(smarts_dict[pattern]) < n_mols:
            if added_to_pattern:
                break
            ss = smarts_search_dict[pattern]
            oechem.OEPrepareSearch(mol, ss)
        
            if ss.SingleMatch(mol):
                oechem.OEAssignFormalCharges(mol)
                net_charge = oechem.OENetCharge(mol)
                if len(smarts_dict[pattern]) == 0:
                    charge_key = f"net_abs_charge_{abs(net_charge)}"
                    smarts_dict[pattern] = {charge_key: []}
                    smarts_dict[pattern][charge_key].append((f"{smiles} Pubchem_CID_{mol_name}", get_tree_fp(mol)))
                    added_to_pattern = True
                    mols_added.append(smiles)
                else:
                    charge_key = f"net_abs_charge_{abs(net_charge)}"
                    if charge_key not in smarts_dict[pattern]:
                        smarts_dict[pattern][charge_key] = []
                    
                    if len(smarts_dict[pattern][charge_key]) < n_mols:
                        fpmol = get_tree_fp(mol)
                        sim_with_rest = [oegraphsim.OETanimoto(fpmol, fp) for _, fp in smarts_dict[pattern][charge_key]]
                        # If tanimoto sim is less than 0.4 then add molecule to the list
                        if all(x < 0.4 for x in sim_with_rest):
                            smarts_dict[pattern][charge_key].append((f"{smiles} Pubchem_CID_{mol_name}", fpmol))
                            added_to_pattern = True
                            mols_added.append(smiles)
        else:
            all_smarts.remove(pattern)

print(all_smarts)
print('###########################\n')
print(set(all_smarts) - set(smarts_dict.keys()))

for key, value in smarts_dict.items():
    smarts_dict[key] = {charge_key: [x[0] for x in molecules] for charge_key, molecules in value.items()}

with open(f'./individual_json_files/{fname[:-7]}_smarts_dict.json', "w") as f:
    json.dump(smarts_dict, f, indent=4)
