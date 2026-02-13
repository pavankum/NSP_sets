# NSP_sets
preprocessing scripts for pubchem search

## File Manifest

| File | Description |
|------|-------------|
| `submit_search.sh` | Shell script to submit PubChem search jobs for the subsets named within|
| `pubchem_NPS_search.sh` | Main script for executing PubChem NPS (nitrogen/phosphorous/sulfur) searches |
| `1.pubchem_oe_pattern_search.py` | Python script for pattern-based molecule searching using OpenEye toolkit |
| `2.merge_all_jsons.py` | Merges multiple JSON result files into a single consolidated file |
| `3.restrict_to_n_heavy_atoms_oeversion.py` | Further filtering of molecules from combined search results using OpenEye |
| `FILTER.OE` | OpenEye filter definition file for molecule selection criteria |
| `nitrogen_summary_updated.md` | Summary document for nitrogen-containing molecule SMARTS |
| `sulfur_summary_updated.md` | Summary document for sulfur-containing molecule SMARTS |
| `phosphorous_summary_updated.md` | Summary document for phosphorous-containing molecule SMARTS |
| `set1-N.json` | JSON results for nitrogen-containing molecules |
| `set1-N-smiles.smi` | SMILES format file for nitrogen subset |
| `set1-P-smiles.smi` | SMILES format file for phosphorous subset |
| `set1-S-smiles.smi` | SMILES format file for sulfur subset |



