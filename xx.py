#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create a Boltz-2-compatible ligand pickle from a minimal mmCIF:
- keeps atom names from _atom_site.label_atom_id
- uses given 3D coords (Cartn_x/y/z)
- infers bonds from coordinates (no CCD bond loop required)
- sets PDB_NAME=<CODE>, symmetries=<hex-pickled permutations>
- pickles with rdkit AllProps to <out_dir>/<CODE>.pkl

Usage:
  python make_boltz_ligand_pkl_from_cif.py RSB.cif /path/to/mols RAX
"""

import sys
import pathlib
import pickle
from typing import List, Dict, Any

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds as RDB
from rdkit.Chem.rdchem import Conformer

# Try Biopython (optional)
try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # type: ignore
    _HAS_BIO = True
except Exception:
    _HAS_BIO = False


def _parse_cif_minimal(cif_path: str) -> Dict[str, List[str]]:
    """
    Minimal loop_ parser for the atom_site table.
    Returns a dict of column_name -> list of values (all strings).
    """
    cols = []
    rows = []
    in_loop = False
    current_cols = []
    data: Dict[str, List[str]] = {}

    with open(cif_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower() == "loop_":
                in_loop = True
                current_cols = []
                continue
            if in_loop and line.startswith("_atom_site."):
                current_cols.append(line)
                continue
            if in_loop and current_cols and not line.startswith("_"):
                # row of data for this loop
                rows.append(line)
                # detect end of loop when we hit a non-data token (handled after loop)
                continue
            # if we hit another token and we were in atom_site loop, stop reading it
            if in_loop and current_cols and (line.startswith("_") or line.lower() == "loop_"):
                break

    if not current_cols or not rows:
        raise ValueError("Could not find an _atom_site loop in the CIF.")

    # split rows by whitespace (CIF is whitespace-separated in your example)
    ncol = len(current_cols)
    table = [[] for _ in range(ncol)]
    for r in rows:
        parts = r.split()
        if len(parts) < ncol:
            # tolerate trailing missing fields by padding
            parts = parts + ["?"] * (ncol - len(parts))
        elif len(parts) > ncol:
            # some CIF writers may put quoted strings with spaces; keep it simple here
            parts = parts[:ncol]
        for i, v in enumerate(parts):
            table[i].append(v)

    for name, col in zip(current_cols, table):
        data[name] = col
    return data


def load_atoms_from_cif(cif_path: str) -> List[Dict[str, Any]]:
    """
    Extract atoms as dicts with: element, name, x, y, z, group_PDB.
    Prefer Biopython; fall back to minimal parser.
    """
    keys = {
        "grp": "_atom_site.group_PDB",
        "sym": "_atom_site.type_symbol",
        "nam": "_atom_site.label_atom_id",
        "x": "_atom_site.Cartn_x",
        "y": "_atom_site.Cartn_y",
        "z": "_atom_site.Cartn_z",
    }

    if _HAS_BIO:
        d = MMCIF2Dict(cif_path)
        # Biopython returns scalars or lists; normalize to lists
        def as_list(v):
            return v if isinstance(v, list) else [v]
        try:
            grp = as_list(d[keys["grp"]])
            sym = as_list(d[keys["sym"]])
            nam = as_list(d[keys["nam"]])
            xs  = [float(x) for x in as_list(d[keys["x"]])]
            ys  = [float(y) for y in as_list(d[keys["y"]])]
            zs  = [float(z) for z in as_list(d[keys["z"]])]
        except KeyError:
            # fallback to minimal
            d = _parse_cif_minimal(cif_path)
            grp = d.get(keys["grp"], [])
            sym = d.get(keys["sym"], [])
            nam = d.get(keys["nam"], [])
            xs  = [float(x) for x in d.get(keys["x"], [])]
            ys  = [float(y) for y in d.get(keys["y"], [])]
            zs  = [float(z) for z in d.get(keys["z"], [])]
    else:
        d = _parse_cif_minimal(cif_path)
        grp = d.get(keys["grp"], [])
        sym = d.get(keys["sym"], [])
        nam = d.get(keys["nam"], [])
        xs  = [float(x) for x in d.get(keys["x"], [])]
        ys  = [float(y) for y in d.get(keys["y"], [])]
        zs  = [float(z) for z in d.get(keys["z"], [])]

    n = min(len(grp), len(sym), len(nam), len(xs), len(ys), len(zs))
    atoms = []
    for i in range(n):
        atoms.append({
            "group_PDB": grp[i],
            "element": sym[i],
            "name": nam[i],
            "x": xs[i],
            "y": ys[i],
            "z": zs[i],
        })
    return atoms


def build_rdkit_from_atoms(atoms: List[Dict[str, Any]], code: str) -> Chem.Mol:
    """Build an RDKit Mol with one conformer from atom list; infer bonds from coords."""
    rw = Chem.RWMol()
    conf = Conformer()
    conf.Set3D(True)

    pt = Chem.GetPeriodicTable()

    idx_map = []
    for i, a in enumerate(atoms):
        if a["group_PDB"] != "HETATM":
            continue
        elem = a["element"]
        try:
            z = int(pt.GetAtomicNumber(elem))
            if z == 0:
                raise ValueError
        except Exception:
            raise ValueError(f"Unknown element in CIF: {elem!r} (row {i+1})")

        rd_a = Chem.Atom(z)
        # store names like CCD does
        name = str(a["name"])
        rd_a.SetProp("name", name)
        rd_a.SetProp("alt_name", name)
        # leaving_atom default to 0
        rd_a.SetProp("leaving_atom", "0")
        ai = rw.AddAtom(rd_a)
        idx_map.append(ai)

        conf.SetAtomPosition(ai, Chem.rdGeometry.Point3D(a["x"], a["y"], a["z"]))

    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)

    # infer bonds from geometry
    RDB.DetermineBonds(mol, charge=0, useHueckel=True, allowChargedFragments=True)

    # sanitize (permissive)
    try:
        Chem.SanitizeMol(mol, catchErrors=True)
    except Exception:
        pass

    # Boltz metadata
    mol.SetProp("PDB_NAME", code)

    # compute symmetries like Boltz script (skip leaving atoms)
    _compute_symmetries(mol)

    return mol


def _compute_symmetries(mol: Chem.Mol):
    """Compute self-permutation symmetries on non-leaving atoms and store hex-pickled list."""
    import pickle as pkl

    m = Chem.RemoveHs(mol, sanitize=False)
    idx_map = {}
    nxt = 0
    for i, a in enumerate(m.GetAtoms()):
        leave = 0
        if a.HasProp("leaving_atom"):
            try:
                leave = int(a.GetProp("leaving_atom"))
            except Exception:
                leave = 0
        if leave:
            continue
        idx_map[i] = nxt
        nxt += 1

    perms = []
    for raw in m.GetSubstructMatches(m, uniquify=False):
        if {raw[i] for i in idx_map} == set(idx_map.keys()):
            perms.append([idx_map[i] for i in raw if i in idx_map])

    mol.SetProp("symmetries", pkl.dumps(perms).hex())


def write_boltz_pickle(mol: Chem.Mol, out_dir: str, code: str):
    # ensure RDKit pickles include all properties
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"{code}.pkl", "wb") as f:
        pickle.dump(mol, f, protocol=pickle.HIGHEST_PROTOCOL)


def main():
    if len(sys.argv) != 4:
        print("Usage: python make_boltz_ligand_pkl_from_cif.py <ligand.cif> <out_mols_dir> <3letter_code>")
        sys.exit(2)

    cif_path, out_dir, code = sys.argv[1], sys.argv[2], sys.argv[3]
    atoms = load_atoms_from_cif(cif_path)
    if not atoms:
        raise RuntimeError("No atoms parsed from CIF.")
    mol = build_rdkit_from_atoms(atoms, code)
    # tag conformer (nice-to-have)
    try:
        c = mol.GetConformer(0)
        c.SetProp("name", "ProvidedCoords")
        c.SetProp("coord_generation", "from_CIF_Cartn")
    except Exception:
        pass
    write_boltz_pickle(mol, out_dir, code)
    print(f"Wrote {pathlib.Path(out_dir)/f'{code}.pkl'}")


if __name__ == "__main__":
    main()
