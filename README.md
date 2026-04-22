# Expression Rescue Pipeline

ProteinMPNN-based analysis to propose CDR single-point mutations that may rescue antibody expression without disrupting antigen binding. Takes one antibody–antigen complex at a time.

**View the example run with full styling** (GitHub strips `<style>` blocks, so pandas dashboards render plain on github.com):

[![View example on nbviewer](https://img.shields.io/badge/nbviewer-view_example_result-F37626?logo=jupyter&logoColor=white)](https://nbviewer.org/github/CSSB-SNU/ab-expression-rescue/blob/main/example/expression_rescue_example_result.ipynb)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CSSB-SNU/ab-expression-rescue/blob/main/example/expression_rescue_example_result.ipynb)

## Overview

Antibodies obtained from *in silico* design, phage-display, or directed evolution can bind their target well yet still express poorly in cell culture. Some CDR residues are preferred for **binding** (stabilized by the antigen contact) while others are preferred only for **expression / solubility** of the antibody alone. Distinguishing the two guides safe rescue mutations: keep the binding residues, change the rest.

This pipeline uses [ProteinMPNN](https://github.com/dauparas/LigandMPNN)'s single-AA conditional logits as a proxy for residue preference. For each input complex we score two states:

- **bound**: the complex as provided
- **unbound**: the same coordinates with the antigen chain(s) removed

Residues where the wild-type amino acid is strongly preferred in the bound state (but not in the unbound state) are **key binding residues**. Residues where another amino acid out-scores the wild-type in the unbound state are **rescue candidates**.

### Key outputs

The **final output is a short list of rescue mutations** (e.g. `V52S, S108L, T100L`) shown in a prominent green panel at the end of Step 6 and stored in `rescue_ranking.csv` under `top_k=True`. Everything else the notebook produces — heatmaps, logit-colored PDBs, per-residue CSVs, inline 3-D views — is **additional information** that helps you understand how those mutations were selected. 

### ⚠ Output quality depends on input structure quality

This pipeline reads per-residue structural context (backbone geometry + local environment) to produce the logits that drive every downstream decision. If the input complex has low resolution, misfolded regions, poor quality, etc., the resulting rescue mutations may be unreliable. Prefer high-resolution crystal/cryo-EM structures, or high-quality predicted structures with energy-minimization whose binding interface has been validated. 

## Directory layout

```
expression_rescue_pipeline/
├── README.md                          # This file
├── expression_rescue.ipynb            # End-to-end pipeline notebook (one structure per run)
├── ProteinMPNN/                       # Minimal scoring slice of LigandMPNN
│   ├── README.md                      # What this is, upstream link, deviations
│   ├── LICENSE                        # MIT
│   ├── requirements.txt               # Python deps for scoring
│   ├── score.py                       # Scoring entry point
│   ├── data_utils.py                  # PDB parsing + featurization
│   ├── model_utils.py                 # ProteinMPNN model definition
│   └── model_params/
│       └── proteinmpnn_v_48_020.pt    # Pre-trained weights (6.4 MB)
└── example/
    ├── HALP_example.pdb                          # Example antibody-antigen complex
    ├── inputs_example.txt                        # Example user-inputs for the notebook
    ├── expression_rescue_example_result.ipynb    # Frozen example run (outputs preserved, read-only reference)
    └── example_output/                           # Artifacts produced by the example run:
        ├── heatmap.png                               #   3-panel logit heatmap
        ├── rescue_ranking.csv                        #   residues ranked by WT unbound logit
        ├── key_binding_residues.csv                  #   Δ = bound − unbound per residue
        ├── HALP_example_unbound.pdb                  #   antigen-removed complex
        ├── HALP_example_logits.pdb                   #   bound PDB colored by bound logit (B-factor)
        └── HALP_example_unbound_logits.pdb           #   unbound PDB colored by unbound logit (B-factor)
```

## Installation

Tested with Python 3.10 and CUDA 12.1.

```bash
conda create -n expression_rescue python=3.10 -y
conda activate expression_rescue
pip install -r ProteinMPNN/requirements.txt
pip install jupyter matplotlib seaborn pandas nbformat py3Dmol "setuptools<81"
```

The model weights ship with the repo under `ProteinMPNN/model_params/proteinmpnn_v_48_020.pt` — no separate download needed.

## Usage

1. Allocate a GPU (The scoring step requires a GPU. CPU-only execution works but it may be slow).
2. Open `expression_rescue.ipynb` in Jupyter.
3. Edit the **User inputs** cell:
   - `bound_pdb_path`: path to the antibody–antigen complex PDB
   - `antigen_chain_ids`: list of chain letters that make up the antigen, e.g. `["A"]`
   - `rescue_residues`: list of `"<chain><resnum>"` strings to evaluate, e.g. `["B30", "B31", ...]`
   - `output_dir`: where results will be written
   - `binding_logit_threshold`: Δ above which Step 5 flags a position as key-binding (default `1.0`)
   - `top_k_rescue`: how many top-ranked positions Step 6 highlights as the mutations to introduce (default `3`)
   - `user_key_binding_residues` *(optional)*: positions you already know are binding-critical from experiments / prior knowledge; leave empty `[]` to rely solely on the pipeline's Δ-based detection
   - `user_key_binding_mode`: `"union"` (default — pipeline ∪ user) or `"user_only"` (trust only user input). Ignored when `user_key_binding_residues` is empty.
4. **Run all cells.** When it finishes, scroll to Step 6: the **green "Final rescue mutations" box** is the final rescue mutations that our pipeline recommends.

The notebook is self-contained: each step generates the inputs for the next, so re-running from the top is safe. 

## Input format

- **PDB file** — standard ATOM records. Chain letters matter: antibody and antigen chains must differ. Multi-model files are not supported (pass the first model only).
- **Residue spec** — `"<chain_letter><resnum>"`, e.g. `"B30"` = residue 30 on chain B. Insertion codes or chain ID with more than 1 character are not supported by this selector; if your PDB uses them, renumber or rechain first.

## Outputs

**You can see major outputs at jupyter notebook output cells in rendered version**, and all raw outputs written to `<output_dir>/`. **The file most users need is `rescue_ranking.csv` (rows with `top_k=True`)**; other outputs are for user's understanding and further analysis if needed.

### Key output (raw output in `<output_dir>/`)

| File | Contents |
|------|----------|
| `rescue_ranking.csv` | All rescue residues ranked by WT unbound logit ascending. **Rows with `top_k=True` are the selected mutations** — the same ones shown in the green Step 6 panel. `is_key_binding=True` rows were skipped. |

### Additional outputs (raw output in `<output_dir>/`)

| File | Contents |
|------|----------|
| `key_binding_residues.csv` | Per-residue Δ table + flags (`is_key_binding`, `user_specified`, `is_final_key_binding`) |
| `heatmap.png` | Three-panel heatmap: bound / unbound / Δ |
| `<stem>_logits.pdb` | Bound PDB with B-factor set to the WT bound logit (open in PyMOL with `spectrum b, red_white_blue`) |
| `<stem>_unbound_logits.pdb` | Unbound PDB with B-factor set to the WT unbound logit |
| `<stem>_unbound.pdb` | Antigen-removed complex used for unbound scoring |
| `bound/<stem>.pt`, `unbound/<stem>.pt` | ProteinMPNN raw scoring output (tensors) |

Step 8 renders both logit-colored structures inline in the notebook via `py3Dmol` using a **red-white-blue** gradient (low logit → red = WT disfavored; high logit → blue = WT preferred). For a higher-fidelity view, open either file in PyMOL:

```
spectrum b, red_white_blue
```

## How to read the outputs

- **`key_binding_residues.csv`** — residues with large positive Δ are binding-critical. Leave these alone.
- **`rescue_ranking.csv`** — every rescue residue listed, sorted by `WT_unbound_logit` ascending. The top rows are the positions where the antibody-alone model is least happy with the current residue (i.e. the WT amino acid is the least preferred at that site). `best_aa` is the amino acid with the highest unbound logit at that position, and `logit_diff = best_unbound_logit − WT_unbound_logit`. The `top_k` column flags the positions selected for mutation.
- **`heatmap.png`** — rows are residues (sorted chain → resnum), columns are the 21-letter alphabet. The red `*` marks the current WT. Dark cells = high logit on the first two panels. The difference panel is diverging red-white-blue centered at 0: **blue** cells next to the WT `*` mark potential binding contacts (antigen favors that AA); **red** cells mark amino acids the antigen disfavors.

### Mutation-selection protocol (used in the paper)

1. Run the pipeline and inspect `rescue_ranking.csv`.
2. **Exclude key-binding residues first.** Step 5 builds the *final* key-binding set from two sources:
   - **Δ-based detection** (always on): positions with Δ = logit_bound − logit_unbound > `binding_logit_threshold` (default 1.0).
   - **User-supplied list** (optional, via `user_key_binding_residues`): positions you already know are binding-critical from experiments, epitope mapping, or prior literature. The pipeline merges them with the Δ-based set according to `user_key_binding_mode`:
     - `"union"` *(default)* — `pipeline ∪ user`
     - `"user_only"` — trust only the user list, ignore Δ detection.

   The final set is stored in `key_binding_residues.csv` under `is_final_key_binding`. All residues in this set are **automatically skipped** when Step 6 picks the top-K. Mutating them would risk the antigen interface, so they are never rescue targets.
3. From the remaining non-binding residues, take the **top 3** by ascending `WT_unbound_logit` (the positions where the WT amino acid is least favored when the antigen is removed — the strongest expression-rescue targets that don't touch the binding interface).
4. At each selected position, introduce the `best_aa` substitution (the amino acid with the highest unbound logit). The notebook prints the mutations in the form `WT{resnum}best_aa` directly below the ranking table, and the same selection is stored in `rescue_ranking.csv` under `top_k=True`.
5. `top_k_rescue` in the user-inputs cell controls the size of the selected set if you want a different cutoff.

## Citations

- *(This repo's companion paper — citation will be added here after publication.)*
- Dauparas, J. et al. *Robust deep learning-based protein sequence design using ProteinMPNN.* **Science** 378, 49–56 (2022).
- Dauparas, J. et al. *Atomic context-conditioned protein sequence design using LigandMPNN.* **Nature Methods** (2025).

## License

- Pipeline code (`expression_rescue.ipynb`, `README.md`, `example/*`): **MIT**, © 2026 Kiheesoo. See `LICENSE`.
- Bundled ProteinMPNN scoring slice (`ProteinMPNN/`): **MIT**, © 2024 Justas Dauparas, redistributed from [LigandMPNN](https://github.com/dauparas/LigandMPNN). See `ProteinMPNN/LICENSE`.
