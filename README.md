# Expression Rescue Pipeline

ProteinMPNN-based analysis to propose CDR single-point mutations that may rescue antibody expression without disrupting antigen binding. Takes one antibody–antigen complex at a time.

## Overview

Antibodies obtained from *in silico* design, phage-display, or directed evolution can bind their target well yet still express poorly in cell culture. Some CDR residues are preferred for **binding** (stabilized by the antigen contact) while others are preferred only for **expression / solubility** of the antibody alone. Distinguishing the two guides safe rescue mutations: keep the binding residues, change the rest.

This pipeline uses [ProteinMPNN](https://github.com/dauparas/LigandMPNN)'s single-AA conditional logits as a proxy for residue preference. For each input complex we score two states:

- **bound**: the complex as provided
- **unbound**: the same coordinates with the antigen chain(s) removed

Residues where the wild-type amino acid is strongly preferred in the bound state (but not in the unbound state) are **key binding residues**. Residues where another amino acid out-scores the wild-type in the unbound state are **rescue candidates**.

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
    ├── bound_example.pdb                         # Example antibody-antigen complex (Adalimumab_variant)
    ├── inputs_example.txt                        # Example user-inputs for the notebook
    ├── expression_rescue_example_result.ipynb    # Frozen example run (outputs preserved, read-only reference)
    └── example_output/                           # Artifacts produced by the example run:
        ├── heatmap.png                               #   3-panel logit heatmap
        ├── rescue_ranking.csv                        #   residues ranked by WT unbound logit
        ├── key_binding_residues.csv                  #   Δ = bound − unbound per residue
        ├── bound_example_unbound.pdb                 #   antigen-removed complex
        ├── bound_example_logits.pdb                  #   bound PDB colored by bound logit (B-factor)
        └── bound_example_unbound_logits.pdb          #   unbound PDB colored by unbound logit (B-factor)
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
4. **Run all cells.**

The notebook is self-contained: each step generates the inputs for the next, so re-running from the top is safe. Intermediate artifacts are cached on disk, so re-running after a tweak to a downstream step does not re-score.

## Input format

- **PDB file** — standard ATOM records. Chain letters matter: antibody and antigen chains must differ. Multi-model files are not supported (pass the first model only).
- **Residue spec** — `"<chain_letter><resnum>"`, e.g. `"B30"` = residue 30 on chain B. Insertion codes are not supported by this selector; if your PDB uses them, renumber first.

## Outputs

All written to `<output_dir>/`:

| File | Contents |
|------|----------|
| `<stem>_unbound.pdb` | Antigen-removed complex used for unbound scoring |
| `bound/<stem>.pt` | ProteinMPNN raw scoring output (bound state) |
| `unbound/<stem>.pt` | ProteinMPNN raw scoring output (unbound state) |
| `key_binding_residues.csv` | All rescue residues with Δ = logit_bound − logit_unbound |
| `rescue_ranking.csv` | All rescue residues ranked by WT unbound logit ascending; `top_k` column flags the positions selected for mutation |
| `heatmap.png` | Three-panel heatmap: bound / unbound / Δ |
| `<stem>_logits.pdb` | Bound PDB with B-factor set to the WT bound logit |
| `<stem>_unbound_logits.pdb` | Unbound PDB with B-factor set to the WT unbound logit |

Step 8 renders both logit-colored structures inline in the notebook via `py3Dmol` using a **red-white-blue** gradient (low logit → red = WT disfavored; high logit → blue = WT preferred) with a shared color scale so the two panels are directly comparable. For a higher-fidelity view, open either file in PyMOL:

```
spectrum b, red_white_blue
```

## How to read the outputs

- **`key_binding_residues.csv`** — residues with large positive Δ are binding-critical. Leave these alone.
- **`rescue_ranking.csv`** — every rescue residue listed, sorted by `WT_unbound_logit` ascending. The top rows are the positions where the antibody-alone model is least happy with the current residue (i.e. the WT amino acid is the least preferred at that site). `best_aa` is the amino acid with the highest unbound logit at that position, and `logit_diff = best_unbound_logit − WT_unbound_logit`. The `top_k` column flags the positions selected for mutation.
- **`heatmap.png`** — rows are residues (sorted chain → resnum), columns are the 21-letter alphabet. The red `*` marks the current WT. Dark cells = high logit on the first two panels. The difference panel is diverging red-white-blue centered at 0: **blue** cells next to the WT `*` mark binding contacts (antigen favors that AA); **red** cells mark amino acids the antigen disfavors.

### Mutation-selection protocol (used in the paper)

1. Run the pipeline and inspect `rescue_ranking.csv`.
2. Take the **top 3** positions by ascending `WT_unbound_logit` (the positions where the WT amino acid is least favored when the antigen is removed — these are the strongest expression-rescue targets).
3. At each of those three positions, introduce the `best_aa` substitution (the amino acid with the highest unbound logit). The notebook prints the three mutations in the form `WT{resnum}best_aa` directly below the ranking table, and the same selection is stored in `rescue_ranking.csv` under `top_k=True`.
4. `top_k_rescue` in the user-inputs cell controls the size of the highlighted set if you want a different cutoff.

## Citations

- *(This repo's companion paper — add citation here after publication.)*
- Dauparas, J. et al. *Robust deep learning-based protein sequence design using ProteinMPNN.* **Science** 378, 49–56 (2022).
- Dauparas, J. et al. *Atomic context-conditioned protein sequence design using LigandMPNN.* **Nature Methods** (2025).

## License

- Pipeline code (`expression_rescue.ipynb`, `README.md`, `example/*`): **MIT**, © 2026 Kiheesoo. See `LICENSE`.
- Bundled ProteinMPNN scoring slice (`ProteinMPNN/`): **MIT**, © 2024 Justas Dauparas, redistributed from [LigandMPNN](https://github.com/dauparas/LigandMPNN). See `ProteinMPNN/LICENSE`.
