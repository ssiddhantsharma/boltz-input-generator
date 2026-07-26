#!/usr/bin/env python3
"""Generate Boltz input files (YAML or FASTA) from PDB structures and FASTA sequences.

Detects proteins, nucleic acids, and ligands, fetches ligand SMILES, and writes
Boltz-ready inputs plus a `boltz_control.csv` manifest.
"""
import argparse
import csv
import os
import re
import warnings
from pathlib import Path

import requests
import yaml

warnings.filterwarnings('ignore')

try:
    from prody import parsePDB
    PRODY_AVAILABLE = True
except ImportError:
    PRODY_AVAILABLE = False

# Amino acids, including the common non-standard residue MSE (selenomethionine).
AMINO_ACIDS = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'MSE': 'M',
}

NUCLEOTIDES = {
    'DA': 'A', 'DT': 'T', 'DG': 'G', 'DC': 'C',
    'A': 'A', 'U': 'U', 'G': 'G', 'C': 'C',
}

NON_LIGAND = {'HOH', 'WAT'}

# User-supplied ligand SMILES, populated via --add_ligand_mapping.
LIGAND_SMILES_MAP = {}


# --- SMILES fetching -------------------------------------------------------

def fetch_smiles_local(ligand_code):
    """Return SMILES from the local mapping if present."""
    return LIGAND_SMILES_MAP.get(ligand_code.upper())


def fetch_smiles_rcsb(ligand_code):
    """Fetch SMILES for a CCD code from the RCSB chemical-component API.

    Prefers the canonical OpenEye descriptor, then any canonical SMILES, then
    any SMILES descriptor at all.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None
        descriptors = r.json().get('pdbx_chem_comp_descriptor', [])
    except (requests.RequestException, ValueError):
        return None

    smiles = [d for d in descriptors if d.get('type', '').startswith('SMILES') and d.get('descriptor')]
    for want in (
        lambda d: d['type'] == 'SMILES_CANONICAL' and d.get('program') == 'OpenEye OEToolkits',
        lambda d: d['type'] == 'SMILES_CANONICAL',
        lambda d: True,
    ):
        for d in smiles:
            if want(d):
                return d['descriptor']
    return None


def fetch_smiles(ligand_code, skip_fetch=False):
    """Resolve a ligand code to SMILES: local mapping first, then RCSB."""
    smiles = fetch_smiles_local(ligand_code)
    if smiles:
        return smiles
    if skip_fetch:
        return None
    return fetch_smiles_rcsb(ligand_code)


def _ligand_entry(chain_id, ligand_code, skip_fetch):
    """Build an entry dict for a ligand, using SMILES if resolvable else the CCD code."""
    smiles = fetch_smiles(ligand_code, skip_fetch)
    return {
        'chain_id': f"{chain_id}_{ligand_code}",
        'sequence': smiles if smiles else ligand_code,
        'entity_type': 'smiles' if smiles else 'ccd',
        'is_nucleic': False,
        'is_ligand': True,
        'ligand_code': ligand_code,
    }


# --- PDB parsing -----------------------------------------------------------

def extract_sequences_prody(pdb_path, skip_fetch=False):
    """Extract protein, nucleic, and ligand entries using ProDy."""
    if not PRODY_AVAILABLE:
        return []
    try:
        pdb = parsePDB(str(pdb_path))
        if pdb is None:
            return []

        entries = []

        protein_sel = pdb.select('protein')
        if protein_sel is not None:
            for chain in set(protein_sel.getChids()):
                chain_sel = protein_sel.select(f'chain {chain} and name CA')
                if chain_sel is not None:
                    seq = [AMINO_ACIDS[r] for r in chain_sel.getResnames() if r in AMINO_ACIDS]
                    if seq:
                        entries.append({
                            'chain_id': chain,
                            'sequence': "".join(seq),
                            'entity_type': 'protein',
                            'is_nucleic': False,
                            'is_ligand': False,
                        })

        for chain in set(pdb.getChids()):
            chain_sel = pdb.select(f"chain {chain}")
            if chain_sel is None:
                continue
            residues = sorted(zip(chain_sel.getResnums(), chain_sel.getResnames()), key=lambda x: x[0])
            nucleic_seq = [NUCLEOTIDES[r] for _, r in residues if r in NUCLEOTIDES]
            if nucleic_seq:
                entries.append({
                    'chain_id': chain,
                    'sequence': "".join(nucleic_seq),
                    'entity_type': 'rna' if 'U' in nucleic_seq else 'dna',
                    'is_nucleic': True,
                    'is_ligand': False,
                })

        ligands = {}
        all_atoms = pdb.select('all')
        if all_atoms is not None:
            for chain in set(all_atoms.getChids()):
                chain_sel = pdb.select(f'chain {chain}')
                if chain_sel is None:
                    continue
                for rname, _ in sorted(set(zip(chain_sel.getResnames(), chain_sel.getResnums())), key=lambda x: x[1]):
                    if rname not in AMINO_ACIDS and rname not in NUCLEOTIDES and rname not in NON_LIGAND:
                        ligands.setdefault(chain, set()).add(rname)

        for chain, resnames in ligands.items():
            for rname in resnames:
                entries.append(_ligand_entry(chain, rname, skip_fetch))

        return entries
    except Exception as e:
        print(f"ProDy extraction failed: {e}")
        return []


def extract_sequences_manual(pdb_path, skip_fetch=False):
    """Parse a PDB file directly (fallback when ProDy is unavailable)."""
    chains_data = {}
    ligands = {}

    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                chain_id = line[21]
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()

                if line.startswith('ATOM'):
                    chains_data.setdefault(chain_id, {}).setdefault(res_num, res_name)
                elif res_name not in AMINO_ACIDS and res_name not in NUCLEOTIDES and res_name not in NON_LIGAND:
                    ligands.setdefault(chain_id, set()).add(res_name)

        entries = []
        for chain_id, residues in chains_data.items():
            sequence = ''
            nucleic_seq = []
            for res_num in sorted(residues):
                res_name = residues[res_num]
                if res_name in AMINO_ACIDS:
                    sequence += AMINO_ACIDS[res_name]
                elif res_name in NUCLEOTIDES:
                    nucleic_seq.append(NUCLEOTIDES[res_name])
                else:
                    sequence += 'X'  # Unknown residue

            if nucleic_seq:
                entries.append({
                    'chain_id': chain_id,
                    'sequence': ''.join(nucleic_seq),
                    'entity_type': 'rna' if 'U' in nucleic_seq else 'dna',
                    'is_nucleic': True,
                    'is_ligand': False,
                })
            elif sequence:
                entries.append({
                    'chain_id': chain_id,
                    'sequence': sequence,
                    'entity_type': 'protein',
                    'is_nucleic': False,
                    'is_ligand': False,
                })

        for chain_id, resnames in ligands.items():
            for rname in resnames:
                entries.append(_ligand_entry(chain_id, rname, skip_fetch))

        return entries
    except Exception as e:
        print(f"Manual parsing failed: {e}")
        return []


def extract_sequences_from_pdb(pdb_path, skip_fetch=False):
    """Extract sequences from a PDB, preferring ProDy and falling back to manual parsing."""
    entries = extract_sequences_prody(pdb_path, skip_fetch)
    if entries:
        return entries
    return extract_sequences_manual(pdb_path, skip_fetch)


# --- FASTA parsing ---------------------------------------------------------

def parse_fasta_header(header):
    """Parse a Boltz-style FASTA header: `>CHAIN|ENTITY[|MSA]`."""
    if header.startswith('>'):
        header = header[1:]

    parts = header.split('|')
    if len(parts) < 2:
        return None

    entity_type = parts[1].strip().lower()
    if entity_type not in ('protein', 'dna', 'rna', 'smiles', 'ccd'):
        return None

    return {
        'chain_id': parts[0].strip(),
        'entity_type': entity_type,
        'msa_path': parts[2].strip() if len(parts) >= 3 else None,
        'is_nucleic': entity_type in ('dna', 'rna'),
        'is_ligand': entity_type in ('smiles', 'ccd'),
    }


def detect_sequence_type(sequence):
    """Guess the entity type of a bare sequence with no Boltz header."""
    sequence = sequence.upper().strip()
    if not sequence:
        return None

    seq_chars = set(sequence)
    if seq_chars.issubset(set('ATCG')):
        return 'dna'
    if seq_chars.issubset(set('AUCG')):
        return 'rna'
    if seq_chars.issubset(set('ACDEFGHIKLMNPQRSTVWY')):
        return 'protein'
    if len(sequence) > 10 and re.match(r'^[A-Za-z0-9@+\-\[\]()=#$%/.\\]+$', sequence):
        return 'smiles'
    if len(sequence) <= 5 and sequence.isalpha():
        return 'ccd'
    return 'protein'


def get_next_chain(used_chains):
    """Return the next unused chain ID (A..Z, then AA..ZZ, then AAA..ZZZ)."""
    from itertools import product
    letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
    for width in (1, 2, 3):
        for combo in product(letters, repeat=width):
            chain = ''.join(combo)
            if chain not in used_chains:
                return chain
    return None


def process_fasta_entry(header, sequence, used_chains, msa_dir=None, force_type=None):
    """Build an entry dict from one FASTA record (Boltz header or auto-detected)."""
    boltz_info = parse_fasta_header(header)
    if boltz_info:
        entry = {
            'chain_id': boltz_info['chain_id'],
            'entity_type': boltz_info['entity_type'],
            'sequence': sequence,
            'is_nucleic': boltz_info['is_nucleic'],
            'is_ligand': boltz_info['is_ligand'],
        }
        if force_type and boltz_info['is_nucleic']:
            entry['entity_type'] = force_type
        if boltz_info['entity_type'] == 'protein':
            if boltz_info['msa_path']:
                entry['msa_path'] = boltz_info['msa_path']
            elif msa_dir:
                msa_file = Path(msa_dir) / f"{boltz_info['chain_id']}.a3m"
                entry['msa_path'] = str(msa_file) if msa_file.exists() else 'empty'
            else:
                entry['msa_path'] = 'empty'
        return entry

    detected_type = detect_sequence_type(sequence)
    if force_type and detected_type in ('dna', 'rna'):
        detected_type = force_type

    next_chain = get_next_chain(used_chains)
    if not next_chain:
        header_clean = header[1:] if header.startswith('>') else header
        print(f"Warning: Unable to assign chain ID for {header_clean}")
        return None

    entry = {
        'chain_id': next_chain,
        'entity_type': detected_type,
        'sequence': sequence,
        'is_nucleic': detected_type in ('dna', 'rna'),
        'is_ligand': detected_type in ('smiles', 'ccd'),
    }
    if detected_type == 'protein':
        msa_path = 'empty'
        if msa_dir:
            msa_file = Path(msa_dir) / f"{next_chain}.a3m"
            if msa_file.exists():
                msa_path = str(msa_file)
        entry['msa_path'] = msa_path
    return entry


def parse_fasta_file(fasta_path, msa_dir=None, force_type=None):
    """Parse a (possibly multi-record) FASTA file into entry dicts."""
    entries = []
    used_chains = set()
    current_header = None
    current_seq = []

    def flush():
        if current_header is None:
            return
        entry = process_fasta_entry(current_header, ''.join(current_seq), used_chains, msa_dir, force_type)
        if entry:
            entries.append(entry)
            used_chains.add(entry['chain_id'])

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                flush()
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
    flush()

    return entries


# --- Output ----------------------------------------------------------------

def create_yaml_content(entries, use_msa_server=False):
    """Build the Boltz YAML dict, merging identical sequences into shared chain IDs."""
    proteins, nucleics, ligands = {}, {}, {}

    for entry in entries:
        etype = entry['entity_type']
        if etype == 'protein':
            key = entry['sequence']
            if key in proteins:
                proteins[key]['chains'].append(entry['chain_id'])
            else:
                proteins[key] = {
                    'chains': [entry['chain_id']],
                    'sequence': entry['sequence'],
                    'msa_path': entry.get('msa_path', 'empty'),
                }
        elif etype in ('dna', 'rna'):
            key = (entry['sequence'], etype)
            bucket = nucleics.setdefault(key, {'chains': [], 'sequence': entry['sequence'], 'type': etype})
            bucket['chains'].append(entry['chain_id'])
        elif etype in ('smiles', 'ccd'):
            key = (entry['sequence'], etype)
            bucket = ligands.setdefault(key, {'chains': [], 'sequence': entry['sequence'], 'type': etype})
            bucket['chains'].append(entry['chain_id'])

    sequences = []
    for data in proteins.values():
        seq_entry = {'protein': {'id': sorted(data['chains']), 'sequence': data['sequence']}}
        if not use_msa_server:
            seq_entry['protein']['msa'] = data['msa_path'] or 'empty'
        sequences.append(seq_entry)

    for data in nucleics.values():
        sequences.append({data['type']: {'id': sorted(data['chains']), 'sequence': data['sequence']}})

    for data in ligands.values():
        ligand = {'id': sorted(data['chains'])}
        ligand['smiles' if data['type'] == 'smiles' else 'ccd'] = data['sequence']
        sequences.append({'ligand': ligand})

    return {'sequences': sequences}


def create_fasta_content(entries, use_msa_server=False):
    """Build Boltz FASTA content from entries."""
    lines = []
    for entry in entries:
        header = f">{entry['chain_id']}|{entry['entity_type']}"
        if entry['entity_type'] == 'protein' and 'msa_path' in entry and not use_msa_server:
            header += f"|{entry['msa_path']}"
        lines.append(header)
        lines.append(entry['sequence'])
    return '\n'.join(lines)


def process_files(input_dir, output_dir, msa_dir=None, force_type=None,
                  output_format='yaml', use_msa_server=False, skip_fetch=False):
    """Process every PDB/FASTA file in `input_dir`, writing Boltz inputs + a control CSV."""
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'boltz_control.csv')

    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'YAML_Path', 'Source_Type', 'Protein_Chains',
                         'Nucleic_Chains', 'Ligand_Chains', 'Total_Chains', 'Ligand_Details'])

        processed_files = 0
        input_path = Path(input_dir)
        all_files = sorted(list(input_path.glob('*.pdb')) + list(input_path.glob('*.fa*')))

        for file_path in all_files:
            try:
                file_id = file_path.stem
                print(f"Processing {file_path.name}...")

                if file_path.suffix == '.pdb':
                    entries = extract_sequences_from_pdb(file_path, skip_fetch)
                    source_type = 'PDB'
                else:
                    entries = parse_fasta_file(file_path, msa_dir, force_type)
                    source_type = 'FASTA'

                if not entries:
                    print(f"  No valid sequences found in {file_path.name}")
                    continue

                # Resolve MSA paths for proteins that don't already carry one.
                for entry in entries:
                    if entry['entity_type'] != 'protein':
                        continue
                    if msa_dir and not use_msa_server:
                        msa_file = Path(msa_dir) / f"{entry['chain_id']}.a3m"
                        entry['msa_path'] = str(msa_file) if msa_file.exists() else 'empty'
                    else:
                        entry.setdefault('msa_path', 'empty')

                if output_format.lower() == 'fasta':
                    output_path = os.path.join(output_dir, f"{file_id}.fasta")
                    with open(output_path, 'w') as f:
                        f.write(create_fasta_content(entries, use_msa_server))
                else:
                    output_path = os.path.join(output_dir, f"{file_id}.yaml")
                    with open(output_path, 'w') as f:
                        yaml.dump(create_yaml_content(entries, use_msa_server),
                                  f, default_flow_style=False, sort_keys=False)

                protein_chains, nucleic_chains, ligand_chains, ligand_details = [], [], [], []
                for entry in entries:
                    etype = entry['entity_type']
                    if etype == 'protein':
                        protein_chains.append(entry['chain_id'])
                    elif etype in ('dna', 'rna'):
                        nucleic_chains.append(f"{entry['chain_id']}:{etype}")
                    elif etype in ('smiles', 'ccd'):
                        ligand_chains.append(f"{entry['chain_id']}:{etype}")
                        if 'ligand_code' in entry:
                            ligand_details.append(f"{entry['chain_id']}:{entry['ligand_code']}")

                total_chains = len(protein_chains) + len(nucleic_chains) + len(ligand_chains)
                writer.writerow([
                    file_id, os.path.abspath(output_path), source_type,
                    ','.join(sorted(protein_chains)), ','.join(sorted(nucleic_chains)),
                    ','.join(sorted(ligand_chains)), total_chains, ','.join(ligand_details),
                ])

                processed_files += 1
                print(f"  {file_id}: {total_chains} chains "
                      f"({len(protein_chains)} protein, {len(nucleic_chains)} nucleic, "
                      f"{len(ligand_chains)} ligand)")
                if ligand_details:
                    print(f"    Ligands: {', '.join(ligand_details)}")

            except Exception as e:
                print(f"  Error processing {file_path.name}: {e}")

        if processed_files:
            print(f"\nSuccessfully processed {processed_files} files")
            print(f"Output directory: {output_dir}")
            print(f"Control file: {csv_path}")
        else:
            print("No files were processed")

    return processed_files


def main():
    parser = argparse.ArgumentParser(description='Generate Boltz input files from PDB and FASTA files')
    parser.add_argument('--input', required=True, help='Input directory containing PDB/FASTA files')
    parser.add_argument('--output', required=True, help='Output directory for Boltz input files')
    parser.add_argument('--msa_dir', help='Directory containing .a3m MSA files')
    parser.add_argument('--format', choices=['yaml', 'fasta'], default='yaml', help='Output format')
    parser.add_argument('--use_msa_server', action='store_true',
                        help='Use the Boltz MSA server instead of local MSA files')
    parser.add_argument('--rna', action='store_true', help='Force nucleic acids to RNA')
    parser.add_argument('--dna', action='store_true', help='Force nucleic acids to DNA')
    parser.add_argument('--add_ligand_mapping', nargs=2, metavar=('CODE', 'SMILES'), action='append',
                        help='Add a custom ligand SMILES mapping (repeatable)')
    parser.add_argument('--skip_ligand_fetch', action='store_true',
                        help='Skip online SMILES fetching (use only local mappings)')

    args = parser.parse_args()

    if args.rna and args.dna:
        parser.error("Cannot specify both --rna and --dna")
    force_type = 'rna' if args.rna else ('dna' if args.dna else None)

    if args.add_ligand_mapping:
        for code, smiles in args.add_ligand_mapping:
            LIGAND_SMILES_MAP[code.upper()] = smiles
            print(f"Added custom mapping: {code.upper()} -> {smiles}")

    print(f"ProDy: {'available' if PRODY_AVAILABLE else 'not installed (using manual PDB parser)'}\n")

    process_files(args.input, args.output, args.msa_dir, force_type,
                  args.format, args.use_msa_server, args.skip_ligand_fetch)


if __name__ == "__main__":
    main()
