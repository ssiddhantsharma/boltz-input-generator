#!/usr/bin/env python3
import os
import sys
import yaml
import json
import csv
from pathlib import Path
import argparse


def load_yaml_file(yaml_path):
    """Load and parse YAML file"""
    try:
        with open(yaml_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading {yaml_path}: {e}")
        return None


def convert_boltz_to_rf3(yaml_data, name):
    """Convert Boltz YAML format to RosettaFold3 JSON format"""
    rf3_entry = {
        "name": name,
        "components": []
    }
    
    # Extract sequences from the YAML
    if 'sequences' in yaml_data:
        for seq_entry in yaml_data['sequences']:
            if 'protein' in seq_entry:
                protein = seq_entry['protein']
                component = {
                    "seq": protein['sequence'],
                    "chain_id": protein['id'][0] if isinstance(protein['id'], list) else protein['id']
                }
                rf3_entry['components'].append(component)
    
    return rf3_entry


def analyze_chains(yaml_data):
    """Analyze chains in the YAML data and return chain information"""
    protein_chains = []
    nucleic_chains = []
    ligand_chains = []
    
    if 'sequences' in yaml_data:
        for seq_entry in yaml_data['sequences']:
            if 'protein' in seq_entry:
                chain_id = seq_entry['protein']['id'][0] if isinstance(seq_entry['protein']['id'], list) else seq_entry['protein']['id']
                protein_chains.append(chain_id)
            elif 'rna' in seq_entry or 'dna' in seq_entry:
                # Handle nucleic acids if present
                chain_key = 'rna' if 'rna' in seq_entry else 'dna'
                chain_id = seq_entry[chain_key]['id'][0] if isinstance(seq_entry[chain_key]['id'], list) else seq_entry[chain_key]['id']
                nucleic_chains.append(chain_id)
            # Add more chain types as needed
    
    return protein_chains, nucleic_chains, ligand_chains


def main():
    parser = argparse.ArgumentParser(description='Convert Boltz YAML files to RosettaFold3 JSON format')
    parser.add_argument('input_dir', help='Directory containing Boltz YAML files')
    parser.add_argument('output_dir', help='Output directory for JSON files and CSV')
    parser.add_argument('--csv-name', default='rf3_control.csv', help='Name of the output CSV file')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all YAML files
    yaml_files = list(input_dir.glob('*.yaml')) + list(input_dir.glob('*.yml'))
    
    if not yaml_files:
        print(f"No YAML files found in {input_dir}")
        sys.exit(1)
    
    # Storage for all converted data
    rf3_json_data = []
    csv_data = []
    
    # Process each YAML file
    for yaml_file in yaml_files:
        print(f"Processing {yaml_file.name}...")
        
        # Load YAML
        yaml_data = load_yaml_file(yaml_file)
        if yaml_data is None:
            continue
        
        # Get base name (without extension)
        base_name = yaml_file.stem
        
        # Convert to RF3 format
        rf3_entry = convert_boltz_to_rf3(yaml_data, base_name)
        rf3_json_data.append(rf3_entry)
        
        # Analyze chains for CSV
        protein_chains, nucleic_chains, ligand_chains = analyze_chains(yaml_data)
        
        # Create JSON output path
        json_output_path = output_dir / f"{base_name}.json"
        
        # Write individual JSON file
        with open(json_output_path, 'w') as f:
            json.dump([rf3_entry], f, indent=2)
        
        # Prepare CSV row
        csv_row = {
            'ID': base_name,
            'JSON_Path': str(json_output_path.absolute()),
            'Source_Type': 'Boltz_YAML',
            'Protein_Chains': ','.join(protein_chains),
            'Nucleic_Chains': ','.join(nucleic_chains),
            'Ligand_Chains': ','.join(ligand_chains),
            'Total_Chains': len(protein_chains) + len(nucleic_chains) + len(ligand_chains)
        }
        csv_data.append(csv_row)
    
    # Write combined JSON file
    combined_json_path = output_dir / 'all_sequences.json'
    with open(combined_json_path, 'w') as f:
        json.dump(rf3_json_data, f, indent=2)
    
    print(f"Combined JSON written to: {combined_json_path}")
    
    # Write CSV file
    csv_path = output_dir / args.csv_name
    if csv_data:
        fieldnames = ['ID', 'JSON_Path', 'Source_Type', 'Protein_Chains', 'Nucleic_Chains', 'Ligand_Chains', 'Total_Chains']
        
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_data)
        
        print(f"Control CSV written to: {csv_path}")
        print(f"Processed {len(csv_data)} YAML files")
        
        # Print summary
        print("\nSummary:")
        for row in csv_data:
            print(f"  {row['ID']}: {row['Total_Chains']} chains ({row['Protein_Chains']})")
    
    else:
        print("No valid YAML files were processed")


if __name__ == "__main__":
    main()
