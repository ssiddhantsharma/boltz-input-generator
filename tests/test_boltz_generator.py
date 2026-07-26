"""Tests for boltz_generator.

All tests run offline: ligand fetching is exercised with skip_fetch=True or a
monkeypatched HTTP layer, so no network access is required.
"""
import csv

import pytest
import yaml

import boltz_generator as bg


# --- Sequence type detection ----------------------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("ATCGATCG", "dna"),
    ("AUCGAUCG", "rna"),
    ("MVTPEGNVSLKQ", "protein"),
    ("CC(=O)Oc1ccccc1C(=O)O", "smiles"),
    ("ZN", "ccd"),
    ("", None),
    # Known ambiguity: a short code made only of amino-acid letters (e.g. the
    # CCD code "ATP" = Ala/Thr/Pro) is read as a protein, since protein is
    # checked before ccd. Use an explicit `>chain|ccd` header to force it.
    ("ATP", "protein"),
])
def test_detect_sequence_type(seq, expected):
    assert bg.detect_sequence_type(seq) == expected


# --- FASTA header parsing --------------------------------------------------

def test_parse_fasta_header_protein():
    info = bg.parse_fasta_header(">A|protein")
    assert info["chain_id"] == "A"
    assert info["entity_type"] == "protein"
    assert info["msa_path"] is None
    assert info["is_nucleic"] is False


def test_parse_fasta_header_with_msa():
    info = bg.parse_fasta_header(">A|protein|/path/to/a.a3m")
    assert info["msa_path"] == "/path/to/a.a3m"


def test_parse_fasta_header_ligand():
    info = bg.parse_fasta_header(">L|ccd")
    assert info["is_ligand"] is True


@pytest.mark.parametrize("header", [">A", ">A|bogus", "no-pipe"])
def test_parse_fasta_header_invalid(header):
    assert bg.parse_fasta_header(header) is None


# --- Chain ID allocation ---------------------------------------------------

def test_get_next_chain_basic():
    assert bg.get_next_chain(set()) == "A"
    assert bg.get_next_chain({"A", "B"}) == "C"


def test_get_next_chain_overflow_to_two_letters():
    used = {chr(i) for i in range(ord("A"), ord("Z") + 1)}
    assert bg.get_next_chain(used) == "AA"


# --- SMILES fetching -------------------------------------------------------

def test_fetch_smiles_local(monkeypatch):
    monkeypatch.setitem(bg.LIGAND_SMILES_MAP, "ZN", "[Zn+2]")
    assert bg.fetch_smiles("zn") == "[Zn+2]"


def test_fetch_smiles_skip(monkeypatch):
    # With skip_fetch and no local mapping, no network call should happen.
    def boom(*a, **k):
        raise AssertionError("network access attempted")
    monkeypatch.setattr(bg.requests, "get", boom)
    assert bg.fetch_smiles("ATP", skip_fetch=True) is None


def test_fetch_smiles_rcsb_prefers_openeye_canonical(monkeypatch):
    class FakeResp:
        status_code = 200
        def json(self):
            return {"pdbx_chem_comp_descriptor": [
                {"type": "SMILES", "program": "ACDLabs", "descriptor": "acd"},
                {"type": "SMILES_CANONICAL", "program": "CACTVS", "descriptor": "cactvs"},
                {"type": "SMILES_CANONICAL", "program": "OpenEye OEToolkits", "descriptor": "openeye"},
            ]}
    monkeypatch.setattr(bg.requests, "get", lambda *a, **k: FakeResp())
    assert bg.fetch_smiles("ATP") == "openeye"


def test_fetch_smiles_rcsb_404(monkeypatch):
    class FakeResp:
        status_code = 404
        def json(self):
            return {}
    monkeypatch.setattr(bg.requests, "get", lambda *a, **k: FakeResp())
    assert bg.fetch_smiles("NOPE") is None


# --- Manual PDB parsing ----------------------------------------------------

PDB_TEXT = """\
ATOM      1  CA  MET A   1      11.1  13.2  10.0  1.00  0.00           C
ATOM      2  CA  GLY A   2      14.0  14.0  10.0  1.00  0.00           C
ATOM      3  CA  ALA B   1      21.0  20.0  10.0  1.00  0.00           C
HETATM    4  O   HOH C 101      30.0  30.0  10.0  1.00  0.00           O
HETATM    5  C1  ATP D 201      40.0  40.0  10.0  1.00  0.00           C
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_extract_manual_proteins_and_ligand(tmp_path):
    pdb = _write(tmp_path, "x.pdb", PDB_TEXT)
    entries = bg.extract_sequences_manual(pdb, skip_fetch=True)
    proteins = {e["chain_id"]: e["sequence"] for e in entries if e["entity_type"] == "protein"}
    assert proteins == {"A": "MG", "B": "A"}
    # Water is dropped; ATP kept as a ligand (CCD fallback since fetch is skipped).
    ligands = [e for e in entries if e["is_ligand"]]
    assert len(ligands) == 1
    assert ligands[0]["ligand_code"] == "ATP"
    assert ligands[0]["entity_type"] == "ccd"
    assert ligands[0]["chain_id"] == "D_ATP"


def test_extract_manual_nucleic(tmp_path):
    text = (
        "ATOM      1  P    DA E   1       0.0   0.0   0.0  1.00  0.00           P\n"
        "ATOM      2  P    DT E   2       0.0   0.0   0.0  1.00  0.00           P\n"
    )
    entries = bg.extract_sequences_manual(_write(tmp_path, "n.pdb", text), skip_fetch=True)
    assert len(entries) == 1
    assert entries[0]["entity_type"] == "dna"
    assert entries[0]["sequence"] == "AT"


# --- FASTA file parsing ----------------------------------------------------

def test_parse_fasta_file_mixed(tmp_path):
    fasta = _write(tmp_path, "s.fasta", ">chainA|protein\nMVTPEG\n>seq2|dna\nATCG\n")
    entries = bg.parse_fasta_file(fasta)
    assert [e["entity_type"] for e in entries] == ["protein", "dna"]
    assert entries[0]["msa_path"] == "empty"


def test_parse_fasta_file_autodetect_assigns_chains(tmp_path):
    fasta = _write(tmp_path, "s.fasta", ">first\nMVTPEG\n>second\nKLQWYF\n")
    entries = bg.parse_fasta_file(fasta)
    assert [e["chain_id"] for e in entries] == ["A", "B"]


def test_process_fasta_entry_force_type():
    entry = bg.process_fasta_entry(">x|dna", "ATCG", set(), force_type="rna")
    assert entry["entity_type"] == "rna"


# --- YAML / FASTA output ---------------------------------------------------

def test_create_yaml_merges_identical_proteins():
    entries = [
        {"chain_id": "A", "entity_type": "protein", "sequence": "MVT", "msa_path": "empty"},
        {"chain_id": "B", "entity_type": "protein", "sequence": "MVT", "msa_path": "empty"},
    ]
    out = bg.create_yaml_content(entries)
    assert len(out["sequences"]) == 1
    assert out["sequences"][0]["protein"]["id"] == ["A", "B"]


def test_create_yaml_ligand_smiles_and_ccd():
    entries = [
        {"chain_id": "A_HEM", "entity_type": "smiles", "sequence": "C1=CC", "is_ligand": True},
        {"chain_id": "B_ATP", "entity_type": "ccd", "sequence": "ATP", "is_ligand": True},
    ]
    seqs = bg.create_yaml_content(entries)["sequences"]
    ligands = [s["ligand"] for s in seqs]
    assert {"id": ["A_HEM"], "smiles": "C1=CC"} in ligands
    assert {"id": ["B_ATP"], "ccd": "ATP"} in ligands


def test_create_yaml_omits_msa_with_server():
    entries = [{"chain_id": "A", "entity_type": "protein", "sequence": "MVT", "msa_path": "empty"}]
    out = bg.create_yaml_content(entries, use_msa_server=True)
    assert "msa" not in out["sequences"][0]["protein"]


def test_create_fasta_content():
    entries = [{"chain_id": "A", "entity_type": "protein", "sequence": "MVT", "msa_path": "empty"}]
    text = bg.create_fasta_content(entries)
    assert text == ">A|protein|empty\nMVT"


# --- End-to-end ------------------------------------------------------------

def test_process_files_end_to_end(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    (in_dir / "mini.pdb").write_text(PDB_TEXT)
    (in_dir / "seqs.fasta").write_text(">chainA|protein\nMVTPEG\n")

    n = bg.process_files(str(in_dir), str(out_dir), skip_fetch=True)
    assert n == 2

    mini = yaml.safe_load((out_dir / "mini.yaml").read_text())
    kinds = [list(s)[0] for s in mini["sequences"]]
    assert kinds.count("protein") == 2 and kinds.count("ligand") == 1

    with open(out_dir / "boltz_control.csv") as f:
        rows = {r["ID"]: r for r in csv.DictReader(f)}
    assert rows["mini"]["Ligand_Details"] == "D_ATP:ATP"
    assert rows["seqs"]["Protein_Chains"] == "chainA"


def test_process_files_fasta_output(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    (in_dir / "p.fasta").write_text(">A|protein\nMVTPEG\n")
    bg.process_files(str(in_dir), str(out_dir), output_format="fasta", skip_fetch=True)
    assert (out_dir / "p.fasta").read_text().startswith(">A|protein")
