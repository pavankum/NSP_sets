# NSP_sets
preprocessing scripts for pubchem search

## File Manifest

### Scripts
| File | Description |
|------|-------------|
| `submit_search.sh` | Shell script to submit PubChem search jobs for the subsets named within |
| `1.pubchem_oe_pattern_search.py` | Python script for pattern-based molecule searching using OpenEye toolkit |
| `2.merge_all_jsons.py` | Merges multiple JSON result files into a single consolidated file |
| `3.restrict_to_n_heavy_atoms_oeversion.py` | Further filtering of molecules from combined search results using OpenEye |
| `FILTER.OE` | OpenEye filter definition file for molecule selection criteria |

### Reference Documents
| File | Description |
|------|-------------|
| `nitrogen_summary_updated.csv` | Summary document for nitrogen-containing molecule SMARTS |
| `sulfur_summary_updated.csv` | Summary document for sulfur-containing molecule SMARTS |
| `phosphorous_summary_updated.csv` | Summary document for phosphorous-containing molecule SMARTS |

### Results
Results and output files are located in the `results/` directory for each set:
- `pubchem_NSP_search.sh` — JSON results contain the merged results of all inidivdual subset files
- `set{n}-{element}.json` — JSON results for {element}-containing molecules
- `set{n}-{element}-smiles.smi` — SMILES format file for {element} subset




