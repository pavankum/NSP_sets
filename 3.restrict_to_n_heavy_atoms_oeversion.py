import json
import pandas as pd
from openeye import oemolprop
from openeye import oechem
from collections import defaultdict
from openff.toolkit.topology import Molecule



n_df = pd.read_csv('nitrogen_summary_updated.csv', delimiter=',')
n_df = n_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

s_df = pd.read_csv('sulfur_summary_updated.csv', delimiter=',')
s_df = s_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

p_df = pd.read_csv('phosphorous_summary_updated.csv', delimiter=',')
p_df = p_df.applymap(lambda x: x.rstrip() if isinstance(x, str) else x)

n_smarts = n_df['SMARTS'].to_list()
s_smarts = s_df['SMARTS'].to_list()
p_smarts = p_df['SMARTS'].to_list()


new_dict1 = defaultdict(list)
new_dict2 = defaultdict(list)

with open('pubchem_NSP_search.json', "r") as f:
    data = json.load(f)

n_heavy = 40
n_mols = 20
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

# ---- HELPER FUNCTIONS ----
def subset_dict(d, key_list):
    return {k: d[k] for k in key_list if k in d}

def dump_smiles(smiles_list, path):
    with open(path, "w") as f:
        f.write("\n".join(smiles_list) + "\n")

def dump_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=4)

mols_added = []
for k, v in data.items():
    nd_val = []
    for smi in v:
        if smi in mols_added:
            continue
        
        if len(nd_val) == n_mols:
            break
        try:
            offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
            mol = offmol.to_openeye()
        except Exception as exc:
            print(f"Skipping SMILES {smi} due to error during parsing/conversion: {exc}")
            continue
        
        if filter_obj(mol):
            nd_val.append(smi)
            mols_added.append(smi)
            
    new_dict1[k] = nd_val[:10]
    new_dict2[k] = nd_val[10:]


set1_N = subset_dict(new_dict1, n_smarts)
set1_S = subset_dict(new_dict1, s_smarts)
set1_P = subset_dict(new_dict1, p_smarts)

set2_N = subset_dict(new_dict2, n_smarts)
set2_S = subset_dict(new_dict2, s_smarts)
set2_P = subset_dict(new_dict2, p_smarts)


# ---- write out ----
dump_json(set1_N, f"set1-N.json")
dump_json(set1_S, f"set1-S.json")
dump_json(set1_P, f"set1-P.json")

dump_json(set2_N, f"set2-N.json")
dump_json(set2_S, f"set2-S.json")
dump_json(set2_P, f"set2-P.json")

# Write SMILES files for each set
for set_dict, set_name in [(set1_N, "set1-N"), (set1_S, "set1-S"), (set1_P, "set1-P"),
                            (set2_N, "set2-N"), (set2_S, "set2-S"), (set2_P, "set2-P")]:
    smiles_list = []
    for v in set_dict.values():
        smiles_list.extend(v)
    dump_smiles(smiles_list, f"{set_name}-smiles.smi")

# ---- WRITE JSON SETS ----
dump_json(new_dict1, f'pubchem_NSP_search_upto_{n_heavy}_hac_{int(n_mols/2)}_mols_set1.json')
dump_json(new_dict2, f'pubchem_NSP_search_upto_{n_heavy}_hac_{int(n_mols/2)}_mols_set2.json')



# ---- WRITE SMILES FILES ----
# --- SET 1 ---
all_smiles_set1 = []
for k, v in new_dict1.items():
    all_smiles_set1.extend(v)

dump_smiles(all_smiles_set1, f"all_smiles_upto_{n_heavy}_and_{int(n_mols/2)}mols_set1.smi")

# --- SET 2 ---
all_smiles_set2 = []
for k, v in new_dict2.items():
    all_smiles_set2.extend(v)

dump_smiles(all_smiles_set2, f"all_smiles_upto_{n_heavy}_and_{int(n_mols/2)}mols_set2.smi")
