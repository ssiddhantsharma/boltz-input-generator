# Boltz Input Generator

Converts PDB files and FASTA sequences into input files for [Boltz](https://github.com/jwohlwend/boltz) structure prediction.

## Installation

```bash
pip install PyYAML
```

Optional: Install PyRosetta for better PDB parsing (academic license required).

## Usage

```bash
python boltz_generator.py --input INPUT_DIR --output OUTPUT_DIR [OPTIONS]
```

### Options

| Flag | Description |
|------|-------------|
| `--input` | Input directory with PDB/FASTA files |
| `--output` | Output directory |
| `--use_msa_server` | Use Boltz MSA server (recommended) |
| `--msa_dir` | Directory with .a3m MSA files |
| `--format` | Output format: yaml (default) or fasta |
| `--rna` | Force nucleic acids to RNA |
| `--dna` | Force nucleic acids to DNA |

## Input Types

- **PDB files**: Multi-chain protein structures
- **Protein FASTA**: Amino acid sequences
- **DNA/RNA FASTA**: Nucleotide sequences
- **SMILES**: Small molecule strings (e.g., `CC(=O)O`)
- **CCD codes**: 3-letter chemical codes (e.g., `ATP`)

## Examples

### Basic Usage
```bash
# Process PDB files
python boltz_generator.py --input ./pdbs --output ./boltz_inputs --use_msa_server

# Mixed input types
python boltz_generator.py --input ./structures --output ./boltz_inputs --use_msa_server
```

### With Local MSAs
```bash
python boltz_generator.py --input ./pdbs --output ./boltz_inputs --msa_dir ./msas
```

## Output

```
output_directory/
├── structure1.yaml
├── structure2.yaml
└── boltz_control.csv
```

### Example YAML
```yaml
sequences:
  - protein:
      id: [A, B]
      sequence: MVTPEGNVSL...
      msa: empty
```

### Example YAML (Protein-Ligand)
```yaml
sequences:
  - protein:
      id: [A]
      sequence: MVTPEGNVSL...
      msa: empty
  - ligand:
      id: [B]
      smiles: 'CC(=O)O'
```

## Running Boltz

### Basic Prediction
```bash
boltz predict structure.yaml --use_msa_server
```

### Enhanced Quality
```bash
boltz predict structure.yaml --use_msa_server --use_potentials
```

### Affinity Prediction (Protein-Ligand Only)
```bash
boltz predict complex.yaml --use_msa_server
```

## SLURM Integration

Use the generated `boltz_control.csv` for batch processing:

```bash
#!/bin/bash
#SBATCH --array=1-N
CSV_FILE="boltz_control.csv"

params=$(get_parameters "$CSV_FILE" "$SLURM_ARRAY_TASK_ID")
eval "$params"
boltz predict $YAML_Path --use_msa_server
```

## Boltz Features

| Feature | FASTA | YAML |
|---------|-------|------|
| Proteins/Nucleic acids | ✓ | ✓ |
| Small molecules | ✓ | ✓ |
| Custom MSA | ✓ | ✓ |
| Modified residues | ✗ | ✓ |
| Covalent bonds | ✗ | ✓ |
| Affinity prediction | ✗ | ✓ |

## Affinity Output

For protein-ligand complexes, Boltz generates affinity predictions:

- `affinity_probability_binary`: Binding probability (0-1, for screening)
- `affinity_pred_value`: Binding affinity as log(IC50) in μM (for optimization)

## Common Issues and Best Practices

### MSA Considerations
- **Use MSA server when possible**: `--use_msa_server` is recommended for most cases
- **Single sequence mode**: Using `msa: empty` hurts model performance but may be necessary for some cases
- **MSA quality**: Better MSAs generally lead to better predictions

### Input Preparation
- **Chain limits**: Very long sequences or many chains may require more GPU memory
- **Sequence validation**: Check that protein sequences use standard amino acids
- **Ligand formats**: Verify SMILES strings are valid; use CCD codes for standard molecules

### Boltz Prediction Options
From the [Boltz documentation](https://github.com/jwohlwend/boltz), key flags to consider:

```bash
# Enhanced quality (slower)
boltz predict input.yaml --use_msa_server --use_potentials

# AlphaFold3-like parameters (much slower, higher quality)
boltz predict input.yaml --use_msa_server --recycling_steps 10 --diffusion_samples 25

# For detailed analysis
boltz predict input.yaml --use_msa_server --write_full_pae --write_full_pde
```


### Understanding Boltz Output
Boltz generates comprehensive output files:

- **Structure files** (`.cif`): 3D coordinates with confidence scores
- **Confidence scores** (`confidence_*.json`): Quality metrics
  - `confidence_score`: Overall quality (0-1)
  - `ptm`/`iptm`: Template modeling scores
  - `complex_plddt`: Average confidence per residue
- **Quality matrices** (`.npz`): PAE, PDE, pLDDT for detailed analysis
- **Affinity predictions** (`affinity_*.json`): Only for protein-ligand complexes
  - `affinity_probability_binary`: Binding probability (0-1, for screening)
  - `affinity_pred_value`: log(IC50) in μM (for optimization)
