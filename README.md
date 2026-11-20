# PFBoost demo

PFBoost is a Streamlit browser for residue-resolved hydrogen-exchange protection factors and HX-derived opening free energies. The full collection contains 414 datasets, including 358 with a recorded DOI or PubMed citation. This demo uses 20 DOI-cited datasets to demonstrate publication filters, experiment comparison, measurements and sequence downloads, and PF-colored 3D structures.

## Demo data

| Item | Count |
| --- | ---: |
| DOI-cited datasets | 20 |
| Distinct DOIs | 14 |
| Parent protein labels | 16 |
| Recorded experimental state labels | 15 |
| State labels in Browse | 6 |
| Bundled PDB files with residue mappings | 16 |
| Residue-level PF rows | 803 |
| Rows with a source-reported PF | 371 |
| Datasets with FASTA available | 20 |

The input measurements are traceable to released Biological Magnetic Resonance Data Bank records. The bundled structures use released PDB coordinates. Cite the original paper or accession listed in the entry. [BMRB usage policy](https://bmrb.io/bmrb/) · [wwPDB usage policy](https://www.wwpdb.org/about/usage-policies)

The read-only snapshot occupies 2,356,773 bytes (2.25 MiB). Dataset folders contain scientific metadata, measurement tables, and coordinates and residue mappings where available. Standardized PFs use the intrinsic exchange-rate parameters of Nguyen et al. (2018). [Nguyen et al.](https://doi.org/10.1007/s13361-018-2021-z)

The PFBoost logo was created by the project author.

## Run locally

Use Python 3.11. From this directory, install the pinned dependencies in a virtual environment and start Streamlit.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address 127.0.0.1
```

Open the localhost URL printed by Streamlit. The 3D viewer uses WebGL and fetches NGL 2.4.0 from jsDelivr when opened. Citations link to public reference sites.

## Streamlit Community Cloud

Set the entrypoint to `app.py`, select the repository’s `main` branch, and choose Python 3.11 under **Advanced settings**. The app uses the bundled data and requires no app secrets. [Streamlit deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy)
