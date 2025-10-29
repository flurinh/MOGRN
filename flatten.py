import os
import shutil
import re
from typing import List, Tuple, Optional

OUTPUT_ROOT = "new_opsins_outputs"
DEST_DIR = "flat_new_opsins_outputs"
CLEAN_DIR = os.path.join(DEST_DIR, "clean")

# Accept files like "Foo_model_0.cif" (case-insensitive .cif)
MODEL_RE = re.compile(r'(?i)(?P<prefix>[^\\/]+)_model_(?P<idx>\d+)\.cif$')
RESNAME_LIG = {"LIG1"}  # tweak if your ligand has alternative names

# -------------------- FS utils --------------------

def to_long_path(path: str) -> str:
    """Convert a path to an extended-length path on Windows."""
    if os.name == "nt":
        abs_path = os.path.abspath(path)
        if not abs_path.startswith("\\\\?\\"):
            return "\\\\?\\" + abs_path
        return abs_path
    return path

def clean_sample_name(folder_name: str) -> str:
    base, ext = os.path.splitext(folder_name)
    return base if ext.lower() == ".yaml" else folder_name

def iter_cif_files(sample_root: str):
    """Yield (abs_path, model_idx) for any '*_model_<n>.cif' under sample_root."""
    for dirpath, _dirnames, filenames in os.walk(sample_root):
        for fn in filenames:
            m = MODEL_RE.search(fn)
            if not m:
                continue
            idx = int(m.group("idx"))
            src = os.path.abspath(os.path.join(dirpath, fn))
            yield src, idx

# -------------------- geometry backends --------------------

_backend = None
def _select_backend():
    """Pick geometry backend once."""
    global _backend
    if _backend is not None:
        return _backend
    try:
        import gemmi  # noqa: F401
        _backend = "gemmi"
        return _backend
    except Exception:
        pass
    try:
        import Bio  # noqa: F401
        from Bio.PDB import MMCIFParser  # noqa: F401
        _backend = "biopython"
        return _backend
    except Exception:
        _backend = None
        return _backend

def _min_distance_gemmi(cif_path: str) -> Optional[float]:
    import gemmi
    st = gemmi.read_structure(cif_path)
    # Use first model (predictions are usually single-model files)
    lig_coords: List[gemmi.Position] = []
    other_coords: List[gemmi.Position] = []
    for model in st:
        for chain in model:
            for res in chain:
                is_lig = (res.name in RESNAME_LIG)
                for atom in res:
                    # skip altloc handling (gemmi picks conformer position)
                    if is_lig:
                        lig_coords.append(atom.pos)
                    else:
                        other_coords.append(atom.pos)
        break  # only first model

    if not lig_coords or not other_coords:
        return None

    # Brute-force min; fast enough. If needed, could use gemmi.NeighborSearch.
    best = float("inf")
    for p in lig_coords:
        for q in other_coords:
            d = p.dist(q)
            if d < best:
                best = d
    return best

def _min_distance_biopython(cif_path: str) -> Optional[float]:
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure_id = os.path.basename(cif_path)
    try:
        st = parser.get_structure(structure_id, cif_path)
    except Exception:
        return None

    lig_coords: List[Tuple[float, float, float]] = []
    other_coords: List[Tuple[float, float, float]] = []

    # Use first model
    model0 = next(iter(st))
    for chain in model0:
        for res in chain:
            # resname is like 'LIG'
            name = res.get_resname()
            is_lig = (name in RESNAME_LIG)
            for atom in res:
                x, y, z = atom.get_coord()
                if is_lig:
                    lig_coords.append((x, y, z))
                else:
                    other_coords.append((x, y, z))

    if not lig_coords or not other_coords:
        return None

    def dist2(a, b):
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return dx*dx + dy*dy + dz*dz

    best2 = float("inf")
    for p in lig_coords:
        for q in other_coords:
            d2 = dist2(p, q)
            if d2 < best2:
                best2 = d2
    return best2 ** 0.5

def min_lig_contact_distance(cif_path: str) -> Optional[float]:
    """Return the minimum distance (Å) between any LIG atom and any non-LIG atom."""
    backend = _select_backend()
    if backend == "gemmi":
        return _min_distance_gemmi(cif_path)
    elif backend == "biopython":
        return _min_distance_biopython(cif_path)
    else:
        raise RuntimeError(
            "No CIF parser available. Please install one of:\n"
            "  pip install gemmi\n"
            "  pip install biopython"
        )

# -------------------- selection logic --------------------

def main():
    os.makedirs(DEST_DIR, exist_ok=True)
    os.makedirs(CLEAN_DIR, exist_ok=True)

    backend = _select_backend()
    if backend is None:
        # Fail early with a clear message
        raise SystemExit(
            "ERROR: Unable to parse CIF. Install one of:\n"
            "  pip install gemmi   (recommended)\n"
            "  pip install biopython"
        )
    print(f"Using CIF backend: {backend}")

    for entry in os.listdir(OUTPUT_ROOT):
        subfolder_path = os.path.join(OUTPUT_ROOT, entry)
        if not os.path.isdir(subfolder_path):
            continue

        sample_name = clean_sample_name(entry)
        print(f"\nSample '{sample_name}': scanning for models under:\n  {subfolder_path}")

        # Collect candidate models for this sample
        candidates = sorted(iter_cif_files(subfolder_path), key=lambda t: t[1])
        if not candidates:
            print("  No CIFs found.")
            continue

        best_model_idx = None
        best_model_path = None
        best_dist = None

        for cif_path, idx in candidates:
            try:
                dmin = min_lig_contact_distance(to_long_path(cif_path))
            except Exception as e:
                print(f"  [model {idx}] ERROR parsing '{cif_path}': {e}")
                continue

            if dmin is None:
                print(f"  [model {idx}] No LIG or no non-LIG atoms — skipped.")
                continue

            print(f"  [model {idx}] min(LIG ↔ non-LIG) = {dmin:.3f} Å  ({cif_path})")

            if (best_dist is None) or (dmin < best_dist):
                best_dist = dmin
                best_model_idx = idx
                best_model_path = cif_path

        if best_model_path is None:
            print("  No eligible models with LIG found.")
            continue

        # Copy chosen model to clean/
        dest_filename = f"{sample_name}_model_{best_model_idx}.cif"
        dest_path = os.path.join(CLEAN_DIR, dest_filename)

        try:
            shutil.copy2(to_long_path(best_model_path), to_long_path(dest_path))
            print(f"  -> Selected model {best_model_idx} (min distance {best_dist:.3f} Å)")
            print(f"  -> Copied to: {dest_path}")
        except Exception as e:
            print(f"  ERROR copying best model '{best_model_path}' -> '{dest_path}': {e}")

if __name__ == "__main__":
    main()
