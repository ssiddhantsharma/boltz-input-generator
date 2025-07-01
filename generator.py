#!/usr/bin/env python3
import os
import yaml
import csv
import argparse
from pathlib import Path
import warnings
import re
warnings.filterwarnings('ignore')

try:
    import pyrosetta
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False

def extract_sequences_from_pdb(pdb_path):
    if not PYROSETTA_AVAILABLE:
        return extract_sequences_manual(pdb_path)
    
    try:
        pose = pyrosetta.pose_from_pdb(str(pdb_path))
        chains_data = []
        
        pdb_info = pose.pdb_info()
        if not pdb_info:
            return []
        
        current_chain = None
        current_sequence = []
        
        for i in range(1, pose.size() + 1):
            residue_chain = pdb_info.chain(i)
            residue_aa = pose.residue(i).name1()
            
            if current_chain is None:
                current_chain = residue_chain
            
            if residue_chain != current_chain:
                if current_sequence:
                    chains_data.append({
                        'chain_id': current_chain,
                        'sequence': ''.join(current_sequence),
                        'entity_type': 'protein',
                        'is_nucleic': False,
                        'is_ligand': False
                    })
                current_chain = residue_chain
                current_sequence = [residue_aa]
            else:
                current_sequence.append(residue_aa)
        
        if current_sequence:
            chains_data.append({
                'chain_id': current_chain,
                'sequence': ''.join(current_sequence),
                'entity_type': 'protein',
                'is_nucleic': False,
                'is_ligand': False
            })
        
        return chains_data
    except Exception:
        return extract_sequences_manual(pdb_path)

def extract_sequences_manual(pdb_path):
    chains_data = {}
    
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    
                    if chain_id not in chains_data:
                        chains_data[chain_id] = {}
                    
                    if res_num not in chains_data[chain_id]:
                        chains_data[chain_id][res_num] = res_name
        
        entries = []
        for chain_id, residues in chains_data.items():
            sequence = ''
            for res_num in sorted(residues.keys()):
                res_name = residues[res_num]
                sequence += aa_map.get(res_name, 'X')
            
            if sequence:
                entries.append({
                    'chain_id': chain_id,
                    'sequence': sequence,
                    'entity_type': 'protein',
                    'is_nucleic': False,
                    'is_ligand': False
                })
        
        return entries
    except Exception:
        return []

def parse_fasta_header(header):
    if header.startswith('>'):
        header = header[1:]
    
    parts = header.split('|')
    if len(parts) >= 2:
        chain_id = parts[0].strip()
        entity_type = parts[1].strip().lower()
        msa_path = parts[2].strip() if len(parts) >= 3 else None
        
        valid_types = ['protein', 'dna', 'rna', 'smiles', 'ccd']
        if entity_type not in valid_types:
            return None
            
        return {
            'chain_id': chain_id,
            'entity_type': entity_type,
            'msa_path': msa_path,
            'is_nucleic': entity_type in ['dna', 'rna'],
            'is_ligand': entity_type in ['smiles', 'ccd']
        }
    return None

def detect_sequence_type(sequence):
    sequence = sequence.upper().strip()
    
    if not sequence:
        return None
    
    dna_bases = set('ATCG')
    rna_bases = set('AUCG')
    protein_aa = set('ACDEFGHIKLMNPQRSTVWY')
    
    seq_chars = set(sequence)
    
    if seq_chars.issubset(dna_bases):
        return 'dna'
    elif seq_chars.issubset(rna_bases):
        return 'rna'
    elif seq_chars.issubset(protein_aa):
        return 'protein'
    elif re.match(r'^[A-Za-z0-9@+\-\[\]()=#$%/.\\]+$', sequence):
        return 'smiles'
    elif len(sequence) <= 5 and sequence.isalpha():
        return 'ccd'
    else:
        return 'protein'

def get_next_chain(used_chains):
    for chain in (chr(i) for i in range(ord('A'), ord('Z')+1)):
        if chain not in used_chains:
            return chain
    for first in (chr(i) for i in range(ord('A'), ord('Z')+1)):
        for second in (chr(i) for i in range(ord('A'), ord('Z')+1)):
            chain = first + second
            if chain not in used_chains:
                return chain
    return None

def parse_fasta_file(fasta_path, msa_dir=None, force_type=None):
    entries = []
    current_header = None
    current_seq = []
    used_chains = set()
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequence = ''.join(current_seq)
                    entry = process_fasta_entry(current_header, sequence, used_chains, msa_dir, force_type)
                    if entry:
                        entries.append(entry)
                        used_chains.add(entry['chain_id'])
                
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
    
    if current_header:
        sequence = ''.join(current_seq)
        entry = process_fasta_entry(current_header, sequence, used_chains, msa_dir, force_type)
        if entry:
            entries.append(entry)
    
    return entries

def process_fasta_entry(header, sequence, used_chains, msa_dir=None, force_type=None):
    boltz_info = parse_fasta_header(header)
    if boltz_info:
        entry = {
            'chain_id': boltz_info['chain_id'],
            'entity_type': boltz_info['entity_type'],
            'sequence': sequence,
            'is_nucleic': boltz_info['is_nucleic'],
            'is_ligand': boltz_info['is_ligand']
        }
        
        if boltz_info['entity_type'] == 'protein':
            if boltz_info['msa_path']:
                entry['msa_path'] = boltz_info['msa_path']
            elif msa_dir:
                msa_file = Path(msa_dir) / f"{boltz_info['chain_id']}.a3m"
                entry['msa_path'] = str(msa_file) if msa_file.exists() else 'empty'
            else:
                entry['msa_path'] = 'empty'
        
        if force_type and boltz_info['is_nucleic']:
            entry['entity_type'] = force_type
            
        return entry
    
    header_clean = header[1:] if header.startswith('>') else header
    
    detected_type = detect_sequence_type(sequence)
    if force_type and detected_type in ['dna', 'rna']:
        detected_type = force_type
    
    next_chain = get_next_chain(used_chains)
    
    entry = {
        'chain_id': next_chain,
        'entity_type': detected_type,
        'sequence': sequence,
        'is_nucleic': detected_type in ['dna', 'rna'],
        'is_ligand': detected_type in ['smiles', 'ccd']
    }
    
    if detected_type == 'protein':
        msa_path = 'empty'
        if msa_dir:
            msa_file = Path(msa_dir) / f"{next_chain}.a3m"
            if msa_file.exists():
                msa_path = str(msa_file)
        entry['msa_path'] = msa_path
    
    return entry

def create_yaml_content(entries, msa_dir=None, use_msa_server=False):
    sequences = []
    protein_sequences = {}
    nucleic_sequences = {}
    ligand_sequences = {}
    
    for entry in entries:
        if entry['entity_type'] == 'protein':
            key = entry['sequence']
            if key not in protein_sequences:
                protein_sequences[key] = {
                    'chains': [entry['chain_id']],
                    'sequence': entry['sequence'],
                    'msa_path': entry.get('msa_path', 'empty')
                }
            else:
                protein_sequences[key]['chains'].append(entry['chain_id'])
        
        elif entry['entity_type'] in ['dna', 'rna']:
            key = (entry['sequence'], entry['entity_type'])
            if key not in nucleic_sequences:
                nucleic_sequences[key] = {
                    'chains': [entry['chain_id']],
                    'sequence': entry['sequence'],
                    'type': entry['entity_type']
                }
            else:
                nucleic_sequences[key]['chains'].append(entry['chain_id'])
        
        elif entry['entity_type'] in ['smiles', 'ccd']:
            key = (entry['sequence'], entry['entity_type'])
            if key not in ligand_sequences:
                ligand_sequences[key] = {
                    'chains': [entry['chain_id']],
                    'sequence': entry['sequence'],
                    'type': entry['entity_type']
                }
            else:
                ligand_sequences[key]['chains'].append(entry['chain_id'])
    
    for protein_data in protein_sequences.values():
        seq_entry = {
            'protein': {
                'id': sorted(protein_data['chains']),
                'sequence': protein_data['sequence']
            }
        }
        
        if not use_msa_server:
            if protein_data['msa_path'] and protein_data['msa_path'] != 'empty':
                seq_entry['protein']['msa'] = protein_data['msa_path']
            else:
                seq_entry['protein']['msa'] = 'empty'
        
        sequences.append(seq_entry)
    
    for nucleic_data in nucleic_sequences.values():
        sequences.append({
            nucleic_data['type']: {
                'id': sorted(nucleic_data['chains']),
                'sequence': nucleic_data['sequence']
            }
        })
    
    for ligand_data in ligand_sequences.values():
        ligand_entry = {
            'ligand': {
                'id': sorted(ligand_data['chains'])
            }
        }
        
        if ligand_data['type'] == 'smiles':
            ligand_entry['ligand']['smiles'] = ligand_data['sequence']
        elif ligand_data['type'] == 'ccd':
            ligand_entry['ligand']['ccd'] = ligand_data['sequence']
        
        sequences.append(ligand_entry)
    
    return {'sequences': sequences}

def create_fasta_content(entries, use_msa_server=False):
    fasta_lines = []
    
    for entry in entries:
        header = f">{entry['chain_id']}|{entry['entity_type']}"
        if entry['entity_type'] == 'protein' and 'msa_path' in entry and not use_msa_server:
            header += f"|{entry['msa_path']}"
        fasta_lines.append(header)
        fasta_lines.append(entry['sequence'])
    
    return '\n'.join(fasta_lines)

def process_files(input_dir, output_dir, msa_dir=None, force_type=None, output_format='yaml', use_msa_server=False):
    os.makedirs(output_dir, exist_ok=True)
    
    if PYROSETTA_AVAILABLE:
        try:
            pyrosetta.init('-ignore_unrecognized_res -ignore_zero_occupancy false -load_PDB_components false -mute all')
        except:
            pass
    
    csv_path = os.path.join(output_dir, 'boltz_control.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'YAML_Path', 'Source_Type', 'Protein_Chains', 'Nucleic_Chains', 'Ligand_Chains', 'Total_Chains'])

        processed_files = 0
        input_path = Path(input_dir)
        
        all_files = list(input_path.glob('*.pdb')) + list(input_path.glob('*.fa*'))
        
        for file_path in all_files:
            try:
                file_id = file_path.stem
                
                if file_path.suffix == '.pdb':
                    entries = extract_sequences_from_pdb(file_path)
                    source_type = 'PDB'
                else:
                    entries = parse_fasta_file(file_path, msa_dir, force_type)
                    source_type = 'FASTA'
                
                if not entries:
                    continue
                
                for entry in entries:
                    if entry['entity_type'] == 'protein':
                        if msa_dir and not use_msa_server:
                            msa_file = Path(msa_dir) / f"{entry['chain_id']}.a3m"
                            entry['msa_path'] = str(msa_file) if msa_file.exists() else 'empty'
                        else:
                            entry['msa_path'] = 'empty'
                
                if output_format.lower() == 'yaml':
                    yaml_content = create_yaml_content(entries, msa_dir, use_msa_server)
                    output_path = os.path.join(output_dir, f"{file_id}.yaml")
                    
                    with open(output_path, 'w') as f:
                        yaml.dump(yaml_content, f, default_flow_style=False, sort_keys=False)
                
                elif output_format.lower() == 'fasta':
                    fasta_content = create_fasta_content(entries, use_msa_server)
                    output_path = os.path.join(output_dir, f"{file_id}.fasta")
                    
                    with open(output_path, 'w') as f:
                        f.write(fasta_content)
                
                protein_chains = []
                nucleic_chains = []
                ligand_chains = []
                
                for entry in entries:
                    if entry['entity_type'] == 'protein':
                        protein_chains.append(entry['chain_id'])
                    elif entry['entity_type'] in ['dna', 'rna']:
                        nucleic_chains.append(f"{entry['chain_id']}:{entry['entity_type']}")
                    elif entry['entity_type'] in ['smiles', 'ccd']:
                        ligand_chains.append(f"{entry['chain_id']}:{entry['entity_type']}")
                
                total_chains = len(protein_chains) + len(nucleic_chains) + len(ligand_chains)
                
                writer.writerow([
                    file_id,
                    os.path.abspath(output_path),
                    source_type,
                    ','.join(sorted(protein_chains)),
                    ','.join(sorted(nucleic_chains)),
                    ','.join(sorted(ligand_chains)),
                    total_chains
                ])
                
                processed_files += 1
                
            except Exception as e:
                pass
        
        print(f"Processed {processed_files} files")

def main():
    parser = argparse.ArgumentParser(description='Generate Boltz input files from PDB and FASTA files')
    parser.add_argument('--input', required=True, help='Input directory containing PDB/FASTA files')
    parser.add_argument('--output', required=True, help='Output directory for Boltz input files')
    parser.add_argument('--msa_dir', help='Directory containing .a3m MSA files')
    parser.add_argument('--format', choices=['yaml', 'fasta'], default='yaml', help='Output format')
    parser.add_argument('--use_msa_server', action='store_true', help='Use MSA server instead of local MSA files')
    parser.add_argument('--rna', action='store_true', help='Force nucleic acids to RNA')
    parser.add_argument('--dna', action='store_true', help='Force nucleic acids to DNA')
    
    args = parser.parse_args()
    
    if args.rna and args.dna:
        parser.error("Cannot specify both --rna and --dna")
    
    force_type = 'rna' if args.rna else ('dna' if args.dna else None)
    
    process_files(args.input, args.output, args.msa_dir, force_type, args.format, args.use_msa_server)

if __name__ == "__main__":
    main()
