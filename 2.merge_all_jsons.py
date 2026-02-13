import json
import os
from collections import defaultdict

# Path to the directory with JSON files
directory = "./individual_json_files"

merged = defaultdict(list)
for filename in os.listdir(directory):
    if filename.endswith(".json"):
        filepath = os.path.join(directory, filename)
        with open(filepath, "r") as f:
            data = json.load(f)
            for k, v in data.items():
                if isinstance(v, list):
                    merged[k].extend(v)
                else:
                    merged[k].append(v)

# Convert back to normal dict
merged = dict(merged)

fname = "pubchem_NSP_search.json"
# Save merged dict to a new JSON file
with open(fname, "w") as f:
    json.dump(merged, f, indent=4)

print(f"Merged JSON written to {fname}")

