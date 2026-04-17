# ProteinMPNN (scoring-only slice)

This folder is a **minimal subset of [LigandMPNN](https://github.com/dauparas/LigandMPNN)** extracted for the expression-rescue pipeline in the parent directory. It keeps only what is needed to *score* an input PDB with the original **ProteinMPNN** backbone model; design/generation and side-chain packing are omitted.

## Contents

| File | Source | Purpose |
|------|--------|---------|
| `score.py` | LigandMPNN | Scoring entry point (single-AA conditional logits) |
| `data_utils.py` | LigandMPNN | PDB parsing + featurization |
| `model_utils.py` | LigandMPNN | `ProteinMPNN` model definition |
| `model_params/proteinmpnn_v_48_020.pt` | LigandMPNN release | Pre-trained weights (6.4 MB) |
| `requirements.txt` | LigandMPNN | Python dependencies |
| `LICENSE` | LigandMPNN | MIT |

## Deviations from upstream

Every deviation is marked inline with a `# NOTE` comment. Two files differ:

**`score.py` — logits export.** Adds mean/std of pre-softmax logits to the saved `.pt` dictionary (upstream exposes only probabilities):

- `out_dict["mean_of_logits"]`
- `out_dict["std_of_logits"]`

These keys are consumed by the expression-rescue notebook (`../expression_rescue.ipynb`) to compare bound vs. unbound logit distributions.

**`model_utils.py` — `single_aa_score` condition flip.** Inside `ProteinMPNN.single_aa_score`, the `use_sequence` branch is inverted relative to upstream (`if use_sequence:` here vs. `if not use_sequence:` upstream). The change makes the flag semantically match `score.py`'s CLI: `--use_sequence 1` now conditions residue *i* on all other native amino acids (leave-one-out), and `--use_sequence 0` scores residue *i* from backbone only.

`data_utils.py` is a verbatim copy of the upstream file.

## Standalone invocation

```bash
python score.py \
    --pdb_path path/to/complex.pdb \
    --out_folder path/to/output_dir \
    --model_type protein_mpnn \
    --checkpoint_protein_mpnn model_params/proteinmpnn_v_48_020.pt \
    --batch_size 10 \
    --number_of_batches 1 \
    --single_aa_score 1 \
    --use_sequence 1
```

The output is a single `<pdb_stem>.pt` file containing per-residue logits/probabilities across the 21-letter alphabet (20 AAs + X for unknown), plus residue-name mappings and the native sequence.

## Citation

If you use this code, please cite the original ProteinMPNN and LigandMPNN works:

- Dauparas et al., *Robust deep learning-based protein sequence design using ProteinMPNN*, Science, 2022.
- Dauparas et al., *Atomic context-conditioned protein sequence design using LigandMPNN*, Nature Methods, 2025.

## License

MIT, inherited from LigandMPNN. See `LICENSE`.
