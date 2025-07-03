#!/usr/bin/env python3
import os
import yaml
import csv
import argparse
from pathlib import Path
import warnings
import re
import requests
warnings.filterwarnings('ignore')

try:
    import pyrosetta
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False

try:
    import pypdb
    PYPDB_AVAILABLE = True
except ImportError:
    PYPDB_AVAILABLE = False

try:
    from pubchempy import get_compounds
    PUBCHEMPY_AVAILABLE = True
except ImportError:
    PUBCHEMPY_AVAILABLE = False

try:
    from prody import parsePDB
    PRODY_AVAILABLE = True
except ImportError:
    PRODY_AVAILABLE = False

# Extended amino acid mapping including non-standard residues
AMINO_ACIDS = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M'  # Selenomethionine
}

# Extended nucleotide mapping
NUCLEOTIDES = {
    'DA':'A','DT':'T','DG':'G','DC':'C',
    'A':'A','U':'U','G':'G','C':'C'
}

# Local ligand SMILES mapping (can be extended by user)
LIGAND_SMILES_MAP = {
    # Add custom ligand mappings here
    # 'HEM': 'C1=C([N-]C(=C1CCC(=O)O)CC2=NC(=CC3=NC(=C(C4=CC(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-])C)CCC(=O)O)C)C(=C3C)CCC(=O)O)C2C)C)C.[Fe+2]',
}

def fetch_smiles_pypdb(ligand_code):
    """Try to fetch SMILES via pypdb"""
    if not PYPDB_AVAILABLE:
        return None
    try:
        chem_desc = pypdb.describe_chemical(ligand_code)
        if "describeHet" in chem_desc and "ligandInfo" in chem_desc["describeHet"] and \
           "ligand" in chem_desc["describeHet"]["ligandInfo"]:
            lig_data = chem_desc["describeHet"]["ligandInfo"]["ligand"]
            if "smiles" in lig_data:
                return lig_data["smiles"]
    except:
        pass
    return None

def fetch_smiles_pdbe(ligand_code):
    """Fetch SMILES from PDBeChem REST API"""
    url = f"https://www.ebi.ac.uk/pdbe/chem/api/describe/{ligand_code}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if ligand_code.upper() in data:
                if len(data[ligand_code.upper()]) > 0 and 'smiles' in data[ligand_code.upper()][0]:
                    return data[ligand_code.upper()][0]['smiles']
    except:
        pass
    return None

def fetch_smiles_pubchem(ligand_code):
    """Try to fetch SMILES from PubChem by searching for the ligand code as a name"""
    if not PUBCHEMPY_AVAILABLE:
        return None
    try:
        compounds = get_compounds(ligand_code, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        pass
    return None

def fetch_smiles_local(ligand_code):
    """Return SMILES from local dictionary if available"""
    return LIGAND_SMILES_MAP.get(ligand_code.upper())

def fetch_smiles(ligand_code):
    """Attempt to fetch SMILES using multiple methods in order"""
    # Try local mapping first (fastest)
    smiles = fetch_smiles_local(ligand_code)
    if smiles:
        return smiles
    
    # Check if we should skip online fetching
    if globals().get('SKIP_LIGAND_FETCH', False):
        return None
    
    # Try PDBeChem (most reliable for PDB ligands)
    smiles = fetch_smiles_pdbe(ligand_code)
    if smiles:
        return smiles
    
    # Try pypdb
    smiles = fetch_smiles_pypdb(ligand_code)
    if smiles:
        return smiles
    
    # Try PubChem (might have false positives)
    smiles = fetch_smiles_pubchem(ligand_code)
    if smiles:
        return smiles
    
    return None

def extract_sequences_prody(pdb_path):
    """Extract sequences using ProDy (more robust chain handling)"""
    if not PRODY_AVAILABLE:
        return []
    
    try:
        pdb = parsePDB(str(pdb_path))
        if pdb is None:
            return []
        
        entries = []
        
        # Extract protein sequences
        protein_sel = pdb.select('protein')
        if protein_sel is not None:
            for chain in set(protein_sel.getChids()):
                chain_sel = protein_sel.select(f'chain {chain} and name CA')
                if chain_sel is not None:
                    seq = []
                    for resname in chain_sel.getResnames():
                        if resname in AMINO_ACIDS:
                            seq.append(AMINO_ACIDS[resname])
                    if seq:
                        entries.append({
                            'chain_id': chain,
                            'sequence': "".join(seq),
                            'entity_type': 'protein',
                            'is_nucleic': False,
                            'is_ligand': False
                        })
        
        # Extract DNA/RNA sequences
        chain_ids = set(pdb.getChids())
        for chain in chain_ids:
            chain_sel = pdb.select(f"chain {chain}")
            if chain_sel is None:
                continue
            
            # Extract residues in order
            resnums = chain_sel.getResnums()
            resnames = chain_sel.getResnames()
            chain_residues = sorted(zip(resnums, resnames), key=lambda x: x[0])
            
            nucleic_seq = []
            for _, rname in chain_residues:
                if rname in NUCLEOTIDES:
                    nucleic_seq.append(NUCLEOTIDES[rname])
            
            if nucleic_seq:
                # Simple heuristic: if contains U, it's RNA, otherwise DNA
                entity_type = 'rna' if 'U' in nucleic_seq else 'dna'
                entries.append({
                    'chain_id': chain,
                    'sequence': "".join(nucleic_seq),
                    'entity_type': entity_type,
                    'is_nucleic': True,
                    'is_ligand': False
                })
        
        # Extract ligands
        ligands = {}
        all_atoms = pdb.select('all')
        if all_atoms is not None:
            for chain in set(all_atoms.getChids()):
                chain_sel = pdb.select(f'chain {chain}')
                if chain_sel is None:
                    continue
                unique_residues = sorted(set(zip(chain_sel.getResnames(), chain_sel.getResnums())), key=lambda x: x[1])
                for rname, rnum in unique_residues:
                    if rname not in AMINO_ACIDS and rname not in NUCLEOTIDES and rname not in ['HOH', 'WAT']:
                        if chain not in ligands:
                            ligands[chain] = set()
                        ligands[chain].add(rname)
        
        # Add ligands to entries
        for chain, resnames in ligands.items():
            for rname in resnames:
                smiles = fetch_smiles(rname)
                if smiles:
                    entries.append({
                        'chain_id': f"{chain}_{rname}",
                        'sequence': smiles,
                        'entity_type': 'smiles',
                        'is_nucleic': False,
                        'is_ligand': True,
                        'ligand_code': rname
                    })
                else:
                    entries.append({
                        'chain_id': f"{chain}_{rname}",
                        'sequence': rname,
                        'entity_type': 'ccd',
                        'is_nucleic': False,
                        'is_ligand': True,
                        'ligand_code': rname
                    })
        
        return entries
        
    except Exception as e:
        print(f"ProDy extraction failed: {e}")
        return []

def extract_sequences_from_pdb(pdb_path):
    """Extract sequences from PDB with fallback methods"""
    # Try ProDy first (most comprehensive)
    entries = extract_sequences_prody(pdb_path)
    if entries:
        return entries
    
    # Fallback to PyRosetta
    if PYROSETTA_AVAILABLE:
        try:
            pose = pyrosetta.pose_from_pdb(str(pdb_path))
            chains_data = []
            
            pdb_info = pose.pdb_info()
            if not pdb_info:
                return extract_sequences_manual(pdb_path)
            
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
            pass
    
    # Final fallback to manual parsing
    return extract_sequences_manual(pdb_path)

def extract_sequences_manual(pdb_path):
    """Manual PDB parsing with better ligand detection"""
    chains_data = {}
    ligands = {}
    
    aa_map = AMINO_ACIDS.copy()
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    
                    if line.startswith('ATOM'):
                        # Regular protein/nucleic acid atoms
                        if chain_id not in chains_data:
                            chains_data[chain_id] = {}
                        
                        if res_num not in chains_data[chain_id]:
                            chains_data[chain_id][res_num] = res_name
                    
                    elif line.startswith('HETATM'):
                        # Potential ligands
                        if res_name not in aa_map and res_name not in NUCLEOTIDES and res_name not in ['HOH', 'WAT']:
                            if chain_id not in ligands:
                                ligands[chain_id] = set()
                            ligands[chain_id].add(res_name)
        
        entries = []
        
        # Process protein/nucleic chains
        for chain_id, residues in chains_data.items():
            sequence = ''
            is_nucleic = False
            nucleic_seq = []
            
            for res_num in sorted(residues.keys()):
                res_name = residues[res_num]
                if res_name in aa_map:
                    sequence += aa_map[res_name]
                elif res_name in NUCLEOTIDES:
                    nucleic_seq.append(NUCLEOTIDES[res_name])
                    is_nucleic = True
                else:
                    sequence += 'X'  # Unknown amino acid
            
            if is_nucleic and nucleic_seq:
                entity_type = 'rna' if 'U' in nucleic_seq else 'dna'
                entries.append({
                    'chain_id': chain_id,
                    'sequence': ''.join(nucleic_seq),
                    'entity_type': entity_type,
                    'is_nucleic': True,
                    'is_ligand': False
                })
            elif sequence:
                entries.append({
                    'chain_id': chain_id,
                    'sequence': sequence,
                    'entity_type': 'protein',
                    'is_nucleic': False,
                    'is_ligand': False
                })
        
        # Process ligands
        for chain_id, resnames in ligands.items():
            for rname in resnames:
                smiles = fetch_smiles(rname)
                if smiles:
                    entries.append({
                        'chain_id': f"{chain_id}_{rname}",
                        'sequence': smiles,
                        'entity_type': 'smiles',
                        'is_nucleic': False,
                        'is_ligand': True,
                        'ligand_code': rname
                    })
                else:
                    entries.append({
                        'chain_id': f"{chain_id}_{rname}",
                        'sequence': rname,
                        'entity_type': 'ccd',
                        'is_nucleic': False,
                        'is_ligand': True,
                        'ligand_code': rname
                    })
        
        return entries
        
    except Exception as e:
        print(f"Manual parsing failed: {e}")
        return []

def parse_fasta_header(header):
    """Parse Boltz-style FASTA headers"""
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
    """Improved sequence type detection"""
    sequence = sequence.upper().strip()
    
    if not sequence:
        return None
    
    dna_bases = set('ATCG')
    rna_bases = set('AUCG')
    protein_aa = set('ACDEFGHIKLMNPQRSTVWY')
    
    seq_chars = set(sequence)
    
    # Check for nucleic acids first
    if seq_chars.issubset(dna_bases):
        return 'dna'
    elif seq_chars.issubset(rna_bases):
        return 'rna'
    elif seq_chars.issubset(protein_aa):
        return 'protein'
    # Check for SMILES patterns
    elif re.match(r'^[A-Za-z0-9@+\-\[\]()=#$%/.\\]+$', sequence) and len(sequence) > 10:
        return 'smiles'
    # Short alphabetic sequences might be CCD codes
    elif len(sequence) <= 5 and sequence.isalpha():
        return 'ccd'
    else:
        # Default to protein for ambiguous cases
        return 'protein'

def get_next_chain(used_chains):
    """Generate next available chain ID"""
    # Single letters first
    for chain in (chr(i) for i in range(ord('A'), ord('Z')+1)):
        if chain not in used_chains:
            return chain
    
    # Two-letter combinations
    for first in (chr(i) for i in range(ord('A'), ord('Z')+1)):
        for second in (chr(i) for i in range(ord('A'), ord('Z')+1)):
            chain = first + second
            if chain not in used_chains:
                return chain
    
    # Three-letter combinations (unlikely to be needed)
    for first in (chr(i) for i in range(ord('A'), ord('Z')+1)):
        for second in (chr(i) for i in range(ord('A'), ord('Z')+1)):
            for third in (chr(i) for i in range(ord('A'), ord('Z')+1)):
                chain = first + second + third
                if chain not in used_chains:
                    return chain
    
    return None

def parse_fasta_file(fasta_path, msa_dir=None, force_type=None):
    """Parse FASTA file with improved header handling"""
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
    """Process individual FASTA entry"""
    # Try to parse Boltz-style header first
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
    
    # Fallback to auto-detection
    header_clean = header[1:] if header.startswith('>') else header
    
    detected_type = detect_sequence_type(sequence)
    if force_type and detected_type in ['dna', 'rna']:
        detected_type = force_type
    
    next_chain = get_next_chain(used_chains)
    if not next_chain:
        print(f"Warning: Unable to assign chain ID for {header_clean}")
        return None
    
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
    """Create YAML content for Boltz"""
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
    
    # Add protein sequences
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
    
    # Add nucleic acid sequences
    for nucleic_data in nucleic_sequences.values():
        sequences.append({
            nucleic_data['type']: {
                'id': sorted(nucleic_data['chains']),
                'sequence': nucleic_data['sequence']
            }
        })
    
    # Add ligand sequences
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
    """Create FASTA content for Boltz"""
    fasta_lines = []
    
    for entry in entries:
        header = f">{entry['chain_id']}|{entry['entity_type']}"
        if entry['entity_type'] == 'protein' and 'msa_path' in entry and not use_msa_server:
            header += f"|{entry['msa_path']}"
        fasta_lines.append(header)
        fasta_lines.append(entry['sequence'])
    
    return '\n'.join(fasta_lines)

def process_files(input_dir, output_dir, msa_dir=None, force_type=None, output_format='yaml', use_msa_server=False):
    """Process all files in input directory"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize PyRosetta if available
    if PYROSETTA_AVAILABLE:
        try:
            pyrosetta.init('-ignore_unrecognized_res -ignore_zero_occupancy false -load_PDB_components false -mute all')
        except:
            pass
    
    csv_path = os.path.join(output_dir, 'boltz_control.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ID', 'YAML_Path', 'Source_Type', 'Protein_Chains', 'Nucleic_Chains', 'Ligand_Chains', 'Total_Chains', 'Ligand_Details'])

        processed_files = 0
        input_path = Path(input_dir)
        
        all_files = list(input_path.glob('*.pdb')) + list(input_path.glob('*.fa*'))
        
        for file_path in all_files:
            try:
                file_id = file_path.stem
                print(f"Processing {file_path.name}...")
                
                if file_path.suffix == '.pdb':
                    entries = extract_sequences_from_pdb(file_path)
                    source_type = 'PDB'
                else:
                    entries = parse_fasta_file(file_path, msa_dir, force_type)
                    source_type = 'FASTA'
                
                if not entries:
                    print(f"  No valid sequences found in {file_path.name}")
                    continue
                
                # Add MSA paths for proteins
                for entry in entries:
                    if entry['entity_type'] == 'protein':
                        if msa_dir and not use_msa_server:
                            msa_file = Path(msa_dir) / f"{entry['chain_id']}.a3m"
                            entry['msa_path'] = str(msa_file) if msa_file.exists() else 'empty'
                        else:
                            entry['msa_path'] = 'empty'
                
                # Create output file
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
                
                # Collect statistics
                protein_chains = []
                nucleic_chains = []
                ligand_chains = []
                ligand_details = []
                
                for entry in entries:
                    if entry['entity_type'] == 'protein':
                        protein_chains.append(entry['chain_id'])
                    elif entry['entity_type'] in ['dna', 'rna']:
                        nucleic_chains.append(f"{entry['chain_id']}:{entry['entity_type']}")
                    elif entry['entity_type'] in ['smiles', 'ccd']:
                        ligand_chains.append(f"{entry['chain_id']}:{entry['entity_type']}")
                        if 'ligand_code' in entry:
                            ligand_details.append(f"{entry['chain_id']}:{entry['ligand_code']}")
                
                total_chains = len(protein_chains) + len(nucleic_chains) + len(ligand_chains)
                
                writer.writerow([
                    file_id,
                    os.path.abspath(output_path),
                    source_type,
                    ','.join(sorted(protein_chains)),
                    ','.join(sorted(nucleic_chains)),
                    ','.join(sorted(ligand_chains)),
                    total_chains,
                    ','.join(ligand_details)
                ])
                
                processed_files += 1
                print(f"  {file_id}: {total_chains} chains ({len(protein_chains)} protein, {len(nucleic_chains)} nucleic, {len(ligand_chains)} ligand)")
                if ligand_details:
                    print(f"    Ligands: {', '.join(ligand_details)}")
                
            except Exception as e:
                print(f"  Error processing {file_path.name}: {e}")
        
        if processed_files > 0:
            print(f"\nSuccessfully processed {processed_files} files")
            print(f"Output directory: {output_dir}")
            print(f"Control file: {csv_path}")
            if use_msa_server:
                print("Ready for Boltz with --use_msa_server flag")
            else:
                print("Ready for Boltz prediction")
        else:
            print("No files were processed")

def main():
    parser = argparse.ArgumentParser(description='Generate Boltz input files from PDB and FASTA files')
    parser.add_argument('--input', required=True, help='Input directory containing PDB/FASTA files')
    parser.add_argument('--output', required=True, help='Output directory for Boltz input files')
    parser.add_argument('--msa_dir', help='Directory containing .a3m MSA files')
    parser.add_argument('--format', choices=['yaml', 'fasta'], default='yaml', help='Output format')
    parser.add_argument('--use_msa_server', action='store_true', help='Use MSA server instead of local MSA files')
    parser.add_argument('--rna', action='store_true', help='Force nucleic acids to RNA')
    parser.add_argument('--dna', action='store_true', help='Force nucleic acids to DNA')
    parser.add_argument('--add_ligand_mapping', nargs=2, metavar=('CODE', 'SMILES'), 
                       action='append', help='Add custom ligand SMILES mapping (can be used multiple times)')
    parser.add_argument('--skip_ligand_fetch', action='store_true', 
                       help='Skip automatic ligand SMILES fetching (use only local mappings)')
    
    args = parser.parse_args()
    
    if args.rna and args.dna:
        parser.error("Cannot specify both --rna and --dna")
    
    force_type = 'rna' if args.rna else ('dna' if args.dna else None)
    
    # Add custom ligand mappings if provided
    if args.add_ligand_mapping:
        for code, smiles in args.add_ligand_mapping:
            LIGAND_SMILES_MAP[code.upper()] = smiles
            print(f"Added custom mapping: {code.upper()} -> {smiles}")
    
    # Show available libraries
    print("Available libraries:")
    print(f"  ProDy: {'✓' if PRODY_AVAILABLE else '✗'}")
    print(f"  PyRosetta: {'✓' if PYROSETTA_AVAILABLE else '✗'}")
    print(f"  pypdb: {'✓' if PYPDB_AVAILABLE else '✗'}")
    print(f"  pubchempy: {'✓' if PUBCHEMPY_AVAILABLE else '✗'}")
    print()
    
    # Set global flag for ligand fetching
    global SKIP_LIGAND_FETCH
    SKIP_LIGAND_FETCH = args.skip_ligand_fetch
    
    process_files(args.input, args.output, args.msa_dir, force_type, args.format, args.use_msa_server)

if __name__ == "__main__":
    main()
