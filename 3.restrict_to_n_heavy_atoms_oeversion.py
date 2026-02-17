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
n_mols_1 = 10
n_mols_2 = 20
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
# Modified structure: new_dict1 stores {pattern_key: {charge_key: [smiles]}}
# and new_dict2 stores {pattern_key: {charge_key: [smiles]}}
new_dict1_nested = defaultdict(lambda: defaultdict(list))
new_dict2_nested = defaultdict(lambda: defaultdict(list))

for pattern_key, charge_dict in data.items():
    set1_count = 0
    set2_count = 0
    
    # SET 1: Prioritize neutral molecules (net_abs_charge_0), then add from other charges
    sorted_charges = sorted(charge_dict.items(), key=lambda x: (x[0] != "net_abs_charge_0", x[0]))
    
    for charge_key, smiles_list in sorted_charges:
        if set1_count >= n_mols_1:
            break
        for smi in smiles_list:
            if smi in mols_added or set1_count >= n_mols_1:
                continue
            try:
                offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
                mol = offmol.to_openeye()
            except Exception as exc:
                print(f"Skipping SMILES {smi} due to error during parsing/conversion: {exc}")
                continue
            
            if filter_obj(mol):
                new_dict1_nested[pattern_key][charge_key].append(smi)
                mols_added.append(smi)
                set1_count += 1
                if set1_count >= n_mols_1:
                    break
    
    # SET 2: Diversity across charge states
    # Iterate up to 2 times to ensure coverage across charge states
    for _ in range(2):
        if set2_count >= n_mols_2:
            break
        for charge_key, smiles_list in sorted_charges:
            if set2_count >= n_mols_2:
                break
            molecules_from_charge = 0
            for smi in smiles_list:
                if smi in mols_added or set2_count >= n_mols_2:
                    break
                try:
                    offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
                    mol = offmol.to_openeye()
                except Exception as exc:
                    print(f"Skipping SMILES {smi} due to error during parsing/conversion: {exc}")
                    continue
                
                if filter_obj(mol):
                    new_dict2_nested[pattern_key][charge_key].append(smi)
                    mols_added.append(smi)
                    set2_count += 1
                    molecules_from_charge += 1
                    # Limit to 1 molecule per charge state per iteration
                    if molecules_from_charge >= 1:
                        break

# Flatten nested dicts for backwards compatibility
new_dict1 = {k: [smi for charge_dict in v.values() for smi in charge_dict] for k, v in new_dict1_nested.items()}
new_dict2 = {k: [smi for charge_dict in v.values() for smi in charge_dict] for k, v in new_dict2_nested.items()}

set1_N = subset_dict(new_dict1, n_smarts)
set1_S = subset_dict(new_dict1, s_smarts)
set1_P = subset_dict(new_dict1, p_smarts)

set2_N = subset_dict(new_dict2, n_smarts)
set2_S = subset_dict(new_dict2, s_smarts)
set2_P = subset_dict(new_dict2, p_smarts)

# ---- write out ----
dump_json(set1_N, f"./results/set1-N.json")
dump_json(set1_S, f"./results/set1-S.json")
dump_json(set1_P, f"./results/set1-P.json")

dump_json(set2_N, f"./results/set2-N.json")
dump_json(set2_S, f"./results/set2-S.json")
dump_json(set2_P, f"./results/set2-P.json")

# Write SMILES files for each set
for set_dict, set_name in [(set1_N, "set1-N"), (set1_S, "set1-S"), (set1_P, "set1-P"),
                            (set2_N, "set2-N"), (set2_S, "set2-S"), (set2_P, "set2-P")]:
    smiles_list = []
    for v in set_dict.values():
        smiles_list.extend(v)
    dump_smiles(smiles_list, f"./results/{set_name}-smiles.smi")

# ---- WRITE JSON SETS ----
dump_json(new_dict1_nested, f'./results/pubchem_NSP_search_upto_{n_heavy}_hac_{int(n_mols_1)}_mols_set1.json')
dump_json(new_dict2_nested, f'./results/pubchem_NSP_search_upto_{n_heavy}_hac_{int(n_mols_2)}_mols_set2.json')

# ---- WRITE SMILES FILES ----
# --- SET 1 ---
all_smiles_set1 = []
for k, v in new_dict1.items():
    all_smiles_set1.extend(v)

dump_smiles(all_smiles_set1, f"./results/all_smiles_upto_{n_heavy}_and_{int(n_mols_1)}mols_set1.smi")

# --- SET 2 ---
all_smiles_set2 = []
for k, v in new_dict2.items():
    all_smiles_set2.extend(v)

dump_smiles(all_smiles_set2, f"./results/all_smiles_upto_{n_heavy}_and_{int(n_mols_2)}mols_set2.smi")
