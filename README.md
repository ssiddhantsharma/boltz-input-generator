# Boltz Input Generator

Converts PDB files and FASTA sequences into input files for [Boltz](https://github.com/jwohlwend/boltz) structure prediction with automatic ligand detection and SMILES fetching.

## Installation

### Required Dependencies
```bash
pip install PyYAML requests
```

### Optional Dependencies (Recommended)
```bash
# For enhanced PDB parsing and ligand detection
pip install prody

# For ligand SMILES fetching
pip install pypdb pubchempy

# For additional PDB parsing (academic license required)
# Install PyRosetta from https://www.pyrosetta.org/
```

## Features

- **Automatic ligand detection** from PDB HETATM records
- **Multi-source SMILES fetching** (PDBeChem, pypdb, PubChem)
- **Robust PDB parsing** with ProDy, PyRosetta, and manual fallbacks
- **Extended residue support** (including selenomethionine MSE)
- **Custom ligand mappings** via command line
- **Better nucleic acid detection** (DNA vs RNA)
- **Comprehensive error handling** and processing reports

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
| `--add_ligand_mapping CODE SMILES` | Add custom ligand SMILES (repeatable) |
| `--skip_ligand_fetch` | Skip automatic SMILES fetching |

## Input Types

- **PDB files**: Multi-chain protein structures with automatic ligand detection
- **Protein FASTA**: Amino acid sequences
- **DNA/RNA FASTA**: Nucleotide sequences (auto-detected or forced)
- **SMILES**: Small molecule strings (e.g., `CC(=O)O`)
- **CCD codes**: 3-letter chemical codes (e.g., `ATP`)

## Examples

### Basic Usage with Automatic Ligand Detection
```bash
# Process PDB files with automatic ligand SMILES fetching
python boltz_generator.py --input ./pdbs --output ./boltz_inputs --use_msa_server

# Mixed input types
python boltz_generator.py --input ./structures --output ./boltz_inputs --use_msa_server
```

### Custom Ligand Mappings
```bash
# Add custom ligand SMILES
python boltz_generator.py --input ./pdbs --output ./boltz_inputs \
  --add_ligand_mapping HEM "C1=C([N-]C(=C1CCC(=O)O)CC2=NC..." \
  --add_ligand_mapping ZN "Zn" \
  --use_msa_server
```

### Local Processing Only
```bash
# Skip online SMILES fetching, use only local mappings
python boltz_generator.py --input ./pdbs --output ./boltz_inputs \
  --skip_ligand_fetch --use_msa_server
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
└── boltz_control.csv  # Enhanced with ligand details
```

### Example YAML (Protein Only)
```yaml
sequences:
  - protein:
      id: [A, B]
      sequence: MVTPEGNVSL...
      msa: empty
```

### Example YAML (Protein-Ligand Complex)
```yaml
sequences:
  - protein:
      id: [A]
      sequence: MVTPEGNVSL...
      msa: empty
  - ligand:
      id: [B_HEM]
      smiles: 'C1=C([N-]C(=C1CCC(=O)O)CC2=NC...'
  - ligand:
      id: [C_ATP]
      ccd: 'ATP'
```

### Enhanced Control File

The `boltz_control.csv` now includes detailed ligand information:

| Column | Description |
|--------|-------------|
| ID | Structure identifier |
| YAML_Path | Path to generated YAML file |
| Source_Type | PDB or FASTA |
| Protein_Chains | Protein chain IDs |
| Nucleic_Chains | DNA/RNA chains with type |
| Ligand_Chains | Ligand chains with type |
| Total_Chains | Total number of chains |
| Ligand_Details | Detailed ligand codes found |

## Library Status Display

The script shows which optional libraries are available:

```
Available libraries:
  ProDy: ✓
  PyRosetta: ✗
  pypdb: ✓
  pubchempy: ✓
```

## Ligand Detection & SMILES Fetching

### Automatic Detection
- Scans PDB HETATM records for non-standard residues
- Excludes water (HOH, WAT) and standard amino acids/nucleotides
- Attempts SMILES fetching from multiple sources in order:
  1. Local custom mappings
  2. PDBeChem API (most reliable for PDB ligands)
  3. pypdb library
  4. PubChem (by name search)

### Custom Mappings
```bash
# Add specific ligand SMILES
--add_ligand_mapping HEM "C1=C([N-]C(=C1CCC(=O)O)..."
--add_ligand_mapping ZN "[Zn+2]"
--add_ligand_mapping MG "[Mg+2]"
```

### Fallback Behavior
- If SMILES not found, uses CCD code format
- Generates unique chain IDs for each ligand instance
- Reports all ligands found in processing output

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

Use the enhanced `boltz_control.csv` for batch processing:

```bash
#!/bin/bash
#SBATCH --array=1-N
CSV_FILE="boltz_control.csv"

params=$(get_parameters "$CSV_FILE" "$SLURM_ARRAY_TASK_ID")
eval "$params"

# Check if ligands are present for affinity prediction
if [[ "$Ligand_Chains" != "" ]]; then
    echo "Processing protein-ligand complex: $ID"
    boltz predict $YAML_Path --use_msa_server
else
    echo "Processing protein-only structure: $ID"
    boltz predict $YAML_Path --use_msa_server
fi
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
| **Auto-ligand detection** | ✗ | ✓ |

## Affinity Output

For protein-ligand complexes, Boltz generates affinity predictions:

- `affinity_probability_binary`: Binding probability (0-1, for screening)
- `affinity_pred_value`: Binding affinity as log(IC50) in μM (for optimization)

## Common Issues and Best Practices

### Ligand Handling
- **Custom mappings**: For proprietary or unusual ligands, use `--add_ligand_mapping`
- **Network issues**: Use `--skip_ligand_fetch` if online SMILES fetching fails
- **Manual review**: Check the ligand details in `boltz_control.csv` to verify correct detection
- **Complex ligands**: Some complex ligands may need manual SMILES validation

### MSA Considerations
- **Use MSA server when possible**: `--use_msa_server` is recommended for most cases
- **Single sequence mode**: Using `msa: empty` hurts model performance but may be necessary for some cases
- **MSA quality**: Better MSAs generally lead to better predictions

### Input Preparation
- **Chain limits**: Very long sequences or many chains may require more GPU memory
- **Sequence validation**: Check that protein sequences use standard amino acids
- **Ligand formats**: Verify SMILES strings are valid; use CCD codes for standard molecules
- **Mixed complexes**: The script handles protein-DNA-ligand complexes automatically

### PDB Processing
- **Multiple parsing methods**: Script tries ProDy → PyRosetta → manual parsing
- **Chain ordering**: Uses residue numbers for proper sequence ordering
- **Non-standard residues**: MSE (selenomethionine) is correctly mapped to methionine
- **Error recovery**: Processing continues even if some files fail

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

## Troubleshooting

### Common Error Messages
- **"ProDy extraction failed"**: Falls back to PyRosetta or manual parsing
- **"No valid sequences found"**: Check PDB file format or FASTA headers
- **"Unable to assign chain ID"**: Too many chains (>676), consider splitting input
- **"No_SMILES_found"**: Ligand SMILES not available online, consider custom mapping
- **Network timeouts**: Use `--skip_ligand_fetch` if online fetching is unreliable
