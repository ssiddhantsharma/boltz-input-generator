# Boltz Input Generator

Turn a folder of PDB structures and FASTA sequences into ready-to-run
[Boltz](https://github.com/jwohlwend/boltz) input files (YAML or FASTA), with
automatic ligand detection and SMILES lookup.

For each input it writes one Boltz input file plus a `boltz_control.csv` manifest
summarising the chains found.

## Install

```bash
pip install pyyaml requests            # core
# optional: more robust PDB parsing
pip install prody
```

Or install the package (adds a `boltz-input-generator` command):

```bash
pip install -e .            # add [prody] and/or [test] for extras
```

## Usage

```bash
python boltz_generator.py --input ./structures --output ./boltz_inputs --use_msa_server
```

| Flag | Description |
|------|-------------|
| `--input` | Input directory with `.pdb`/`.fasta` files (required) |
| `--output` | Output directory (required) |
| `--use_msa_server` | Let Boltz generate MSAs at predict time (recommended) |
| `--msa_dir` | Directory of local `.a3m` files, matched by chain ID |
| `--format` | `yaml` (default) or `fasta` |
| `--rna` / `--dna` | Force detected nucleic acids to RNA or DNA |
| `--add_ligand_mapping CODE SMILES` | Supply SMILES for a ligand code (repeatable) |
| `--skip_ligand_fetch` | Don't fetch SMILES online; use local mappings only |

### Examples

```bash
# Custom SMILES for a specific ligand, skip everything else online
python boltz_generator.py --input ./pdbs --output ./out \
  --add_ligand_mapping ZN "[Zn+2]" --skip_ligand_fetch

# Local MSAs instead of the MSA server
python boltz_generator.py --input ./pdbs --output ./out --msa_dir ./msas
```

## Inputs

- **PDB files** — multi-chain proteins, nucleic acids, and HETATM ligands are
  detected automatically (water is ignored; `MSE` maps to methionine).
- **FASTA files** — Boltz-style headers `>CHAIN|ENTITY[|MSA]` are respected
  (`ENTITY` = `protein`/`dna`/`rna`/`smiles`/`ccd`). Headerless records get an
  auto-assigned chain ID and a best-effort type guess.

## Output

```
output/
├── structure1.yaml        # or .fasta with --format fasta
├── structure2.yaml
└── boltz_control.csv       # ID, paths, per-type chain lists, ligand details
```

Example YAML for a protein–ligand complex:

```yaml
sequences:
  - protein:
      id: [A]
      sequence: MVTPEGNVSL...
      msa: empty
  - ligand:
      id: [B_HEM]
      smiles: 'C1=C([N-]...'
  - ligand:
      id: [C_ATP]
      ccd: 'ATP'
```

Then run Boltz on the generated files, e.g.:

```bash
boltz predict structure1.yaml --use_msa_server
```

## Ligand SMILES

For each non-standard HETATM residue the code is resolved to SMILES in order:

1. Local mappings from `--add_ligand_mapping`
2. The [RCSB chemical-component API](https://data.rcsb.org)

If nothing is found, the raw CCD code is emitted (`ccd:` in the YAML) so Boltz
can still resolve it. Use `--add_ligand_mapping` for proprietary or unusual
ligands, or `--skip_ligand_fetch` to stay fully offline.

## Notes

- `msa: empty` runs a protein in single-sequence mode; `--use_msa_server` (or a
  real MSA via `--msa_dir`) gives better predictions.
- Type detection for headerless FASTA is heuristic: a short code made only of
  amino-acid letters (e.g. `ATP`) is read as a protein. Use an explicit
  `>chain|ccd` header to force a ligand.

## Tests

```bash
pip install pytest
pytest
```

Tests run fully offline (SMILES fetching is mocked).

## License

MIT — see [LICENSE](LICENSE).
