import json
import os
from collections import defaultdict

# Path to the directory with JSON files
directory = "./individual_json_files"

merged = defaultdict(lambda: defaultdict(list))
for filename in os.listdir(directory):
    if filename.endswith(".json"):
        filepath = os.path.join(directory, filename)
        with open(filepath, "r") as f:
            data = json.load(f)
            for pattern_key, charge_dict in data.items():
                for charge_key, mol_list in charge_dict.items():
                    merged[pattern_key][charge_key].extend(mol_list)

# Convert back to normal dict
merged = {k: dict(v) for k, v in merged.items()}

fname = "pubchem_NSP_search.json"
# Save merged dict to a new JSON file
with open(fname, "w") as f:
    json.dump(merged, f, indent=4)

print(f"Merged JSON written to {fname}")

