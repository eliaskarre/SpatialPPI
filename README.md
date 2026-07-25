# SpatialPPI

SpatialPPI builds a **spatially constrained 3D representation** of a human protein–protein interaction (PPI) network by combining interaction edges with subcellular localization annotations. The repository contains both:

1. a reusable Python CLI pipeline for coordinate generation, and
2. analysis/plotting notebooks-style scripts and exported figures/tables used during project exploration.

<img width="700" alt="Screenshot 2026-07-25 180736" src="https://github.com/user-attachments/assets/5e37ba49-5131-4900-ab44-276fae24a05b" />


## What this project does

At a high level, the project takes three inputs:

- a PPI edge list,
- a protein localization table,
- a JSON cell-layout configuration.

Then it:

1. builds a global PPI graph,
2. splits the network into per-compartment subgraphs,
3. computes 3D coordinates for each compartment with a constrained force-directed layout,
4. merges all positions into one export table (`positions_per_location.tsv`),
5. optionally visualizes the full cell in interactive 3D (Plotly).

The implementation supports multi-localized proteins (a protein can appear in multiple compartments), per-compartment geometry and repulsion tuning, and deterministic runs via seed control.

---

## What was done in this repository

### 1) Data preparation and mapping

- Scripts in `data/` map gene/protein identifiers to UniProt and normalize localization tables.
- Intermediate and failure-tracking TSVs are stored (e.g., mapped/unmapped IDs and cleaned localization files).
- A consensus PPI edgelist and derived helper tables are included for downstream layout and plotting.

### 2) Constrained cell layout engine

A modular Python implementation in `scripts/` was built with these components:

- **CLI (`scripts/cli.py`)**: argument parsing, filtering (`--only`, `--exclude`), reproducibility controls (`--seed`, `--iterations`), and optional plotting.
- **Config system (`scripts/config.py` + `scripts/config.json`)**: typed location specs and runtime construction of geometric constraints.
- **Constraint and optimizer module (`scripts/constraints.py`)**:
  - shape constraints: sphere, shell, ellipsoid, ellipsoid shell, cylinder,
  - soft-wall penalty forces,
  - custom 3D constrained Fruchterman–Reingold-style spring layout.
- **Model layer (`scripts/model.py`)**: `CellModel` + `Location` abstractions, per-location subgraph extraction, coordinate assignment, and Plotly visualization.
- **I/O helpers (`scripts/io_helpers.py`)**: loading graph/localization inputs and writing a DataDiVR-friendly TSV with RGB and alpha attributes.

### 3) Cell geometry design for 13 compartments

<img width="150" height="150" alt="Screenshot 2026-01-23 172917" src="https://github.com/user-attachments/assets/87293a4c-877c-49f2-bc32-25b5033d4247" />
<img width="90" height="150" alt="Screenshot 2026-01-234" src="https://github.com/user-attachments/assets/2a141a51-ff16-4b13-bcd3-d56a5754bdf9" />
<img width="90" height="150" alt="Screenshot 2026-01-232172917" src="https://github.com/user-attachments/assets/99b77a2b-5ed1-42c0-b2c1-7a86d04a7e68" />

The checked-in default config models 13 subcellular compartments, including cytosol, ER, Golgi, centrosome, actin/microtubule/intermediate filaments, plasma membrane, nucleoplasm, nuclear membrane, mitochondria, primary cilium, and nucleoli.

Each compartment has independent:

- geometric boundary parameters,
- center translation in global cell coordinates,
- minimum degree filter,
- repulsion strength.

### 4) Evaluation plots and exports

The `plots/` folder shows analysis and visualization outputs produced during development, including:

- PPI degree-distribution histograms,
- network snapshots,
- location frequency plots,
- per-compartment PPI renderings,
- summary CSV tables.

### 5) Documentation

- `doc/USER_DOCUMENTATION.md` explains practical execution and input/output expectations.
- `doc/DEVELOPER_DOCUMENTATION.md` explains architecture, data flow, and layout mechanics.

---

## Repository structure

- `scripts/` – main pipeline code
- `data/` – raw/processed data and mapping scripts
- `plots/` – analysis scripts and generated figures/tables
- `doc/` – user and developer documentation
- `labels/` – project label metadata

---

## Quick start

### Requirements

Python 3.12+ with:

- `numpy`
- `pandas`
- `networkx`
- `plotly`

### Run

```bash
python scripts/cli.py \
  --ppi data/consensus_ppi_bioplex_biogrid_intact_huri_edgelist.tsv \
  --localizations data/location_with_uniprot.tsv \
  --config scripts/config.json \
  --outdir output/
```

Optional flags:

- `--only mitochondria,nucleoplasm`
- `--exclude cytosol`
- `--seed 123`
- `--iterations 500`
- `--plot --plot-title "SpatialPPI cell" --no-edges`

---

## Main output

`positions_per_location.tsv` with columns:

- `Node ID`
- `x`, `y`, `z`
- `r`, `g`, `b`, `a`
- `compartment`

This format is designed for downstream visualization workflows (including DataDiVR-style graph ingestion).

---

## Notes and current scope

- This repository already contains both the implementation and several generated artifacts (figures, tables, intermediate TSVs).
- Some scripts are exploratory and path-sensitive (expecting certain working directories/files).
- The core maintained entry point for reproducible runs is `scripts/cli.py` with `scripts/config.json`.

---

## Reference

Pirch, S., Müller, F., Iofinova, E. et al. *The VRNetzer platform enables interactive network analysis in Virtual Reality*. Nat Commun 12, 2432 (2021). https://doi.org/10.1038/s41467-021-22570-w
