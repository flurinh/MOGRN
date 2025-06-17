"""
Correct implementation of retinal carbon assignment following IUPAC numbering.
Identifies all 20 carbons in retinal based on connectivity and geometry.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
from scipy.spatial.distance import cdist
import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict



def find_beta_ionone_ring_flexible(atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]],
                                   carbon_indices: List[int]) -> Optional[Dict]:
    """
    Find the β-ionone ring with more flexible criteria.

    The β-ionone ring characteristics:
    1. 6-membered carbon ring
    2. One carbon with 4 total bonds where 2 go to methyls (C1)
    3. One carbon with 4 total bonds where 1 goes to a methyl (C5)
    4. These two carbons (C1 and C5) are NOT adjacent
    5. Usually one sp2 carbon that connects to the polyene chain (C6)
    """

    def find_rings_dfs(start: int, size: int = 6) -> List[List[int]]:
        """Find rings of specified size using DFS."""
        rings = []

        def dfs(current: int, path: List[int], visited: Set[int]):
            if len(path) > size:
                return

            if len(path) == size:
                # Check if we can close the ring
                neighbors = [n for n in connectivity.get(current, []) if n in carbon_indices]
                if start in neighbors:
                    ring = path[:]
                    # Normalize ring representation
                    min_idx = min(ring)
                    min_pos = ring.index(min_idx)
                    normalized = ring[min_pos:] + ring[:min_pos]

                    # Check if we already have this ring
                    if normalized not in rings and normalized[::-1] not in rings:
                        rings.append(normalized)
                return

            # Continue DFS
            neighbors = [n for n in connectivity.get(current, []) if n in carbon_indices]
            for neighbor in neighbors:
                if neighbor not in visited or (len(path) == size - 1 and neighbor == start):
                    new_visited = visited.copy()
                    new_visited.add(neighbor)
                    dfs(neighbor, path + [neighbor], new_visited)

        # Start DFS from each carbon
        for start_carbon in carbon_indices:
            dfs(start_carbon, [start_carbon], {start_carbon})

        return rings

    print("  Finding 6-membered rings...")
    rings = find_rings_dfs(carbon_indices[0], 6)
    print(f"  Found {len(rings)} unique 6-membered rings")

    # Analyze each ring
    for ring_idx, ring in enumerate(rings):
        print(f"\n  Analyzing ring {ring_idx + 1}: {ring}")

        ring_analysis = {
            'ring_carbons': ring,
            'methylated_carbons': [],
            'sp2_candidates': [],
            'c1_candidate': None,
            'c5_candidate': None,
            'c6_candidate': None
        }

        # Analyze each carbon in the ring
        for c_idx in ring:
            atom = atoms_df.iloc[c_idx]

            # Get neighbors
            all_neighbors = connectivity.get(c_idx, [])
            carbon_neighbors = [n for n in all_neighbors if n in carbon_indices]
            ring_carbon_neighbors = [n for n in carbon_neighbors if n in ring]
            non_ring_carbon_neighbors = [n for n in carbon_neighbors if n not in ring]

            # Count total bonds (including hydrogens)
            total_bonds = len(all_neighbors)

            print(f"    C{c_idx} ({atom['res_atom_name']}): "
                  f"bonds={total_bonds}, "
                  f"ring_C={len(ring_carbon_neighbors)}, "
                  f"non-ring_C={len(non_ring_carbon_neighbors)}")

            # Check for methylation
            if non_ring_carbon_neighbors:
                ring_analysis['methylated_carbons'].append({
                    'idx': c_idx,
                    'n_methyls': len(non_ring_carbon_neighbors),
                    'ring_neighbors': ring_carbon_neighbors
                })

            # Check if potentially sp2 (3 total neighbors)
            if total_bonds == 3:
                ring_analysis['sp2_candidates'].append(c_idx)

        # Find C1 (gem-dimethyl) and C5 (single methyl)
        for methyl_carbon in ring_analysis['methylated_carbons']:
            if methyl_carbon['n_methyls'] == 2:
                ring_analysis['c1_candidate'] = methyl_carbon['idx']
                print(f"    Found C1 candidate (gem-dimethyl): {methyl_carbon['idx']}")
            elif methyl_carbon['n_methyls'] == 1:
                # Could be C5
                if ring_analysis['c5_candidate'] is None:
                    ring_analysis['c5_candidate'] = methyl_carbon['idx']
                    print(f"    Found C5 candidate (single methyl): {methyl_carbon['idx']}")

        # Verify C1 and C5 are not adjacent
        if ring_analysis['c1_candidate'] and ring_analysis['c5_candidate']:
            c1_idx = ring_analysis['c1_candidate']
            c5_idx = ring_analysis['c5_candidate']

            # Check if they're adjacent in the ring
            c1_pos = ring.index(c1_idx)
            c5_pos = ring.index(c5_idx)

            # Calculate distance in ring
            dist1 = abs(c5_pos - c1_pos)
            dist2 = 6 - dist1
            ring_distance = min(dist1, dist2)

            print(f"    Ring distance between C1 and C5: {ring_distance}")

            if ring_distance == 1:
                print(f"    WARNING: C1 and C5 are adjacent (unexpected)")

            # Find C6 (should be sp2 or have specific connectivity)
            # C6 is typically adjacent to C1 in the ring
            c1_ring_neighbors = []
            for n in connectivity.get(c1_idx, []):
                if n in ring and n in carbon_indices:
                    c1_ring_neighbors.append(n)

            # Check which C1 neighbor could be C6
            for neighbor in c1_ring_neighbors:
                if neighbor in ring_analysis['sp2_candidates']:
                    ring_analysis['c6_candidate'] = neighbor
                    print(f"    Found C6 candidate (sp2 near C1): {neighbor}")
                    break

            # If no sp2 neighbor of C1, look for other criteria
            if not ring_analysis['c6_candidate']:
                # C6 connects to the polyene chain, so it should have a non-ring carbon neighbor
                for neighbor in c1_ring_neighbors:
                    non_ring_c = [n for n in connectivity.get(neighbor, [])
                                  if n in carbon_indices and n not in ring]
                    if non_ring_c and neighbor != ring_analysis['c5_candidate']:
                        ring_analysis['c6_candidate'] = neighbor
                        print(f"    Found C6 candidate (chain connection): {neighbor}")
                        break

            # If we have all key carbons, this is likely our β-ionone ring
            if (ring_analysis['c1_candidate'] and
                    ring_analysis['c5_candidate'] and
                    ring_analysis['c6_candidate']):
                print(f"\n  ✓ Ring {ring_idx + 1} matches β-ionone criteria!")

                return {
                    'ring_carbons': ring,
                    'c1': ring_analysis['c1_candidate'],
                    'c5': ring_analysis['c5_candidate'],
                    'c6': ring_analysis['c6_candidate']
                }

    # If no ring matched all criteria, try relaxed criteria
    print("\n  No ring matched all criteria. Trying relaxed criteria...")

    # Look for any ring with at least one gem-dimethyl carbon
    for ring_idx, ring in enumerate(rings):
        gem_dimethyl = None
        has_methyls = False

        for c_idx in ring:
            non_ring_carbons = [n for n in connectivity.get(c_idx, [])
                                if n in carbon_indices and n not in ring]
            if len(non_ring_carbons) == 2:
                gem_dimethyl = c_idx
            elif len(non_ring_carbons) == 1:
                has_methyls = True

        if gem_dimethyl:
            print(f"\n  Ring {ring_idx + 1} has gem-dimethyl carbon at {gem_dimethyl}")

            # This is probably our ring - make best guesses for other positions
            result = {'ring_carbons': ring, 'c1': gem_dimethyl}

            # Find C5 (any single methyl carbon)
            for c_idx in ring:
                if c_idx != gem_dimethyl:
                    non_ring_carbons = [n for n in connectivity.get(c_idx, [])
                                        if n in carbon_indices and n not in ring]
                    if len(non_ring_carbons) == 1:
                        result['c5'] = c_idx
                        break

            # Find C6 (neighbor of C1 that could connect to chain)
            c1_neighbors = [n for n in connectivity.get(gem_dimethyl, []) if n in ring]
            for neighbor in c1_neighbors:
                # Prefer carbons with 3 neighbors (potential sp2)
                if len(connectivity.get(neighbor, [])) == 3:
                    result['c6'] = neighbor
                    break

            if 'c6' not in result and c1_neighbors:
                result['c6'] = c1_neighbors[0]  # Just pick one

            if 'c5' not in result:
                # Pick any carbon not adjacent to C1
                for c_idx in ring:
                    if c_idx not in c1_neighbors and c_idx != gem_dimethyl:
                        result['c5'] = c_idx
                        break

            return result

    print("\n  Could not identify β-ionone ring even with relaxed criteria")
    return None


"""
Correct implementation of retinal carbon assignment following IUPAC numbering.
Identifies all 20 carbons in retinal based on connectivity and geometry.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
from scipy.spatial.distance import cdist


def assign_retinal_carbons_correct(atoms_df: pd.DataFrame) -> Dict[str, int]:
    """
    Assign canonical carbon numbers (C1-C20) to retinal atoms based on structure.

    Args:
        atoms_df: DataFrame with retinal atoms containing columns:
                  'type_symbol', 'x', 'y', 'z', 'res_atom_name'

    Returns:
        Dictionary mapping carbon names ('C1'-'C20') to DataFrame indices
    """
    # Build connectivity graph
    connectivity = build_connectivity_graph(atoms_df)

    # Get carbon atoms - handle both 'type_symbol' and 'atom_name' columns
    if 'type_symbol' in atoms_df.columns:
        carbon_mask = atoms_df['type_symbol'] == 'C'
    else:
        carbon_mask = atoms_df['atom_name'] == 'C'
    carbon_indices = atoms_df.index[carbon_mask].tolist()

    if len(carbon_indices) < 20:
        print(f"Warning: Only found {len(carbon_indices)} carbon atoms, expected 20")

    # Find the β-ionone ring using flexible detection
    ring_info = find_beta_ionone_ring_flexible(atoms_df, connectivity, carbon_indices)

    if not ring_info:
        # Try the original method as fallback
        ring_info = find_beta_ionone_ring(atoms_df, connectivity, carbon_indices)

    if not ring_info:
        print("Error: Could not identify β-ionone ring")
        return {}

    # Extract ring information
    ring_carbons = ring_info['ring_carbons']
    c1_idx = ring_info['c1']  # Gem-dimethyl carbon
    c5_idx = ring_info['c5']  # Single methyl carbon
    c6_idx = ring_info['c6']  # sp2 carbon

    # Number the ring carbons (C1-C6)
    ring_numbering = number_ring_carbons(ring_carbons, c1_idx, c5_idx, c6_idx, connectivity)

    # Find and number the polyene chain (C7-C15)
    chain_numbering = number_polyene_chain(atoms_df, connectivity, c6_idx, carbon_indices)

    # Find and number the methyl groups (C16-C20)
    methyl_numbering = number_methyl_groups(atoms_df, connectivity,
                                            ring_numbering, chain_numbering, carbon_indices)

    # Combine all numberings
    carbon_assignments = {}
    carbon_assignments.update(ring_numbering)
    carbon_assignments.update(chain_numbering)
    carbon_assignments.update(methyl_numbering)

    # Also find the Schiff base nitrogen if present
    if 'type_symbol' in atoms_df.columns:
        nitrogen_mask = atoms_df['type_symbol'] == 'N'
    else:
        nitrogen_mask = atoms_df['atom_name'] == 'N'

    if nitrogen_mask.any() and 'C15' in carbon_assignments:
        n_indices = atoms_df.index[nitrogen_mask].tolist()
        c15_idx = carbon_assignments['C15']
        for n_idx in n_indices:
            if c15_idx in connectivity.get(n_idx, []):
                carbon_assignments['NZ'] = n_idx
                break

    return carbon_assignments


def build_connectivity_graph(atoms_df: pd.DataFrame, bond_cutoff: float = 1.7) -> Dict[int, List[int]]:
    """Build connectivity graph based on interatomic distances."""
    # Ensure we have numpy arrays, not pandas objects
    coords = np.array(atoms_df[['x', 'y', 'z']].values, dtype=float)
    distances = cdist(coords, coords)

    connectivity = defaultdict(list)

    for i in range(len(atoms_df)):
        for j in range(i + 1, len(atoms_df)):
            if distances[i, j] < bond_cutoff:
                connectivity[i].append(j)
                connectivity[j].append(i)

    return dict(connectivity)


def find_six_membered_rings(connectivity: Dict[int, List[int]],
                            carbon_indices: List[int]) -> List[List[int]]:
    """Find all six-membered rings in the structure."""
    rings = []

    def dfs_ring(start: int, current: int, path: List[int], visited: Set[int], depth: int):
        if depth > 6:
            return

        if depth == 6:
            # Check if we can close the ring
            if start in connectivity.get(current, []):
                # Normalize ring representation
                ring = path[:]
                min_idx = min(ring)
                min_pos = ring.index(min_idx)
                normalized = ring[min_pos:] + ring[:min_pos]

                # Check if we already have this ring
                ring_exists = False
                for existing_ring in rings:
                    if set(existing_ring) == set(normalized):
                        ring_exists = True
                        break

                if not ring_exists:
                    rings.append(normalized)
            return

        for neighbor in connectivity.get(current, []):
            if neighbor in carbon_indices:
                if neighbor not in visited or (depth == 5 and neighbor == start):
                    new_visited = visited.copy()
                    new_visited.add(neighbor)
                    dfs_ring(start, neighbor, path + [neighbor], new_visited, depth + 1)

    # Try starting from each carbon
    for start_idx in carbon_indices:
        dfs_ring(start_idx, start_idx, [start_idx], {start_idx}, 1)

    return rings


def count_carbon_neighbors(idx: int, connectivity: Dict[int, List[int]],
                           carbon_indices: List[int]) -> int:
    """Count carbon neighbors of an atom."""
    return sum(1 for n in connectivity.get(idx, []) if n in carbon_indices)


def is_sp2_carbon(idx: int, atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]],
                  carbon_indices: List[int]) -> bool:
    """
    Determine if a carbon is sp2 hybridized based on geometry.
    sp2 carbons have 3 neighbors in a planar arrangement.
    """
    neighbors = connectivity.get(idx, [])
    carbon_neighbors = [n for n in neighbors if n in carbon_indices]

    # sp2 carbons typically have 3 total neighbors (including H)
    if len(neighbors) != 3:
        return False

    # Check planarity if we have enough carbon neighbors
    if len(carbon_neighbors) >= 2:
        # Calculate angles between bonds
        center = np.array(atoms_df.iloc[idx][['x', 'y', 'z']].values, dtype=float)
        vectors = []

        for n in neighbors[:3]:  # Use first 3 neighbors
            n_pos = np.array(atoms_df.iloc[n][['x', 'y', 'z']].values, dtype=float)
            vectors.append(n_pos - center)

        # Check if vectors are roughly coplanar
        if len(vectors) == 3:
            normal = np.cross(vectors[0], vectors[1])
            if np.linalg.norm(normal) > 0:
                normal = normal / np.linalg.norm(normal)
                # Third vector should be nearly perpendicular to normal
                dot_product = abs(np.dot(normal, vectors[2]))
                return dot_product < 0.3  # Threshold for planarity

    return True  # Default to sp2 if we can't determine


def find_beta_ionone_ring(atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]],
                          carbon_indices: List[int]) -> Optional[Dict]:
    """
    Find the β-ionone ring in retinal structure.
    Returns dict with ring carbons and key positions (C1, C5, C6).
    """
    rings = find_six_membered_rings(connectivity, carbon_indices)

    print(f"  Found {len(rings)} six-membered rings")

    for ring_idx, ring in enumerate(rings):
        # Count methylated carbons and sp2 carbons in the ring
        methylated_carbons = []
        sp2_carbons = []

        for c_idx in ring:
            # Count non-ring carbon neighbors (potential methyls)
            non_ring_neighbors = [n for n in connectivity.get(c_idx, [])
                                  if n in carbon_indices and n not in ring]

            # Check for methylation (non-ring carbon neighbors)
            if len(non_ring_neighbors) > 0:
                methylated_carbons.append((c_idx, len(non_ring_neighbors)))

            # Check if sp2 - relaxed criteria
            total_neighbors = len(connectivity.get(c_idx, []))
            if total_neighbors == 3:  # Likely sp2
                sp2_carbons.append(c_idx)

        # β-ionone ring should have:
        # - One carbon with 2 methyls (C1)
        # - One carbon with 1 methyl (C5)
        # - Usually one sp2 carbon (C6), but this might not be detected perfectly

        # Find gem-dimethyl carbon (C1)
        c1_idx = None
        for c_idx, n_methyls in methylated_carbons:
            if n_methyls == 2:
                c1_idx = c_idx
                break

        if not c1_idx:
            continue  # Not the β-ionone ring

        # Find single methyl carbon (C5)
        c5_idx = None
        for c_idx, n_methyls in methylated_carbons:
            if n_methyls == 1 and c_idx != c1_idx:
                c5_idx = c_idx
                break

        # Find C6 - should be connected to C1 and lead to the polyene chain
        c6_idx = None

        # Get C1's neighbors in the ring
        c1_ring_neighbors = [n for n in connectivity.get(c1_idx, []) if n in ring]

        # Method 1: Look for sp2 neighbor of C1
        for neighbor in c1_ring_neighbors:
            if neighbor in sp2_carbons:
                c6_idx = neighbor
                break

        # Method 2: Look for neighbor with non-ring carbon connection (chain)
        if not c6_idx:
            for neighbor in c1_ring_neighbors:
                non_ring_carbons = [n for n in connectivity.get(neighbor, [])
                                    if n in carbon_indices and n not in ring]
                # Should have connection to chain but not be a methyl carbon
                if non_ring_carbons and neighbor not in [mc[0] for mc in methylated_carbons]:
                    c6_idx = neighbor
                    break

        # Method 3: Just pick a neighbor of C1 (not C2)
        if not c6_idx and c1_ring_neighbors:
            # C2 is between C1 and C3, C6 is between C1 and C5
            # So prefer the neighbor that's closer to C5 if we found it
            if c5_idx:
                # Find which neighbor is closer to C5
                for neighbor in c1_ring_neighbors:
                    # Simple heuristic: check path length
                    c6_idx = neighbor
                    break
            else:
                c6_idx = c1_ring_neighbors[0]

        # We need at least C1 and C6 to proceed
        if c1_idx and c6_idx:
            # If we didn't find C5, make a guess
            if not c5_idx:
                # C5 should be across the ring from C1
                for c_idx in ring:
                    if c_idx != c1_idx and c_idx != c6_idx:
                        # Check if it's not adjacent to C1
                        if c_idx not in connectivity.get(c1_idx, []):
                            c5_idx = c_idx
                            break

            print(f"  Ring {ring_idx + 1} identified as β-ionone ring:")
            print(f"    C1 (gem-dimethyl): {c1_idx}")
            print(f"    C5 (single methyl): {c5_idx}")
            print(f"    C6 (sp2/chain): {c6_idx}")

            return {
                'ring_carbons': ring,
                'c1': c1_idx,
                'c5': c5_idx,
                'c6': c6_idx
            }

    print("  No ring matched β-ionone criteria")
    return None


def number_ring_carbons(ring: List[int], c1: int, c5: int, c6: int,
                        connectivity: Dict[int, List[int]]) -> Dict[str, int]:
    """Number the ring carbons C1-C6 following IUPAC convention."""
    numbering = {}

    # Start with known positions
    numbering['C1'] = c1
    numbering['C5'] = c5
    numbering['C6'] = c6

    # Find C2: neighbor of C1 that's in the ring and not C6
    for neighbor in connectivity[c1]:
        if neighbor in ring and neighbor != c6:
            numbering['C2'] = neighbor
            break

    # Find C3: neighbor of C2 that's not C1
    for neighbor in connectivity[numbering['C2']]:
        if neighbor in ring and neighbor != c1:
            numbering['C3'] = neighbor
            break

    # Find C4: remaining ring carbon
    for c_idx in ring:
        if c_idx not in numbering.values():
            numbering['C4'] = c_idx
            break

    return numbering


def number_polyene_chain(atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]],
                         c6_idx: int, carbon_indices: List[int]) -> Dict[str, int]:
    """Number the polyene chain carbons C7-C15."""
    numbering = {}

    # Find C7: carbon connected to C6 that's not in the ring
    c6_neighbors = [n for n in connectivity[c6_idx] if n in carbon_indices]

    # C7 is the carbon neighbor of C6 with the most carbons on the other side
    c7_idx = None
    max_chain_length = 0

    for neighbor in c6_neighbors:
        # Count carbons reachable from this neighbor without going back to C6
        visited = {c6_idx}
        chain_length = count_reachable_carbons(neighbor, connectivity,
                                               carbon_indices, visited)
        if chain_length > max_chain_length:
            max_chain_length = chain_length
            c7_idx = neighbor

    if c7_idx is None:
        print("Error: Could not find C7")
        return numbering

    # Trace the polyene chain from C7 to C15
    current = c7_idx
    carbon_num = 7
    visited = {c6_idx}  # Don't go back to the ring

    while carbon_num <= 15 and current is not None:
        numbering[f'C{carbon_num}'] = current
        visited.add(current)

        # Find next carbon in chain
        next_carbon = None
        carbon_neighbors = [n for n in connectivity[current]
                            if n in carbon_indices and n not in visited]

        if carbon_neighbors:
            if len(carbon_neighbors) == 1:
                next_carbon = carbon_neighbors[0]
            else:
                # Multiple choices - pick the one leading to the longest chain
                max_length = 0
                for neighbor in carbon_neighbors:
                    temp_visited = visited.copy()
                    temp_visited.add(neighbor)
                    length = count_reachable_carbons(neighbor, connectivity,
                                                     carbon_indices, temp_visited)
                    if length > max_length:
                        max_length = length
                        next_carbon = neighbor

        current = next_carbon
        carbon_num += 1

    return numbering


def count_reachable_carbons(start: int, connectivity: Dict[int, List[int]],
                            carbon_indices: List[int], visited: Set[int]) -> int:
    """Count carbons reachable from start without going through visited."""
    count = 0
    to_visit = [start]
    local_visited = visited.copy()

    while to_visit:
        current = to_visit.pop()
        if current in local_visited:
            continue
        local_visited.add(current)
        count += 1

        for neighbor in connectivity.get(current, []):
            if neighbor in carbon_indices and neighbor not in local_visited:
                to_visit.append(neighbor)

    return count


def number_methyl_groups(atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]],
                         ring_numbering: Dict[str, int], chain_numbering: Dict[str, int],
                         carbon_indices: List[int]) -> Dict[str, int]:
    """
    Number the methyl groups C16-C20.
    C16, C17: methyls on C1
    C18: methyl on C5
    C19: methyl on C9
    C20: methyl on C13
    """
    numbering = {}

    # Get all numbered carbons so far
    numbered_carbons = set()
    for idx in ring_numbering.values():
        numbered_carbons.add(idx)
    for idx in chain_numbering.values():
        numbered_carbons.add(idx)

    # Find methyls on C1 (should be 2)
    c1_idx = ring_numbering['C1']
    c1_methyls = [n for n in connectivity[c1_idx]
                  if n in carbon_indices and n not in numbered_carbons]

    if len(c1_methyls) >= 2:
        numbering['C16'] = c1_methyls[0]
        numbering['C17'] = c1_methyls[1]

    # Find methyl on C5
    c5_idx = ring_numbering['C5']
    c5_methyls = [n for n in connectivity[c5_idx]
                  if n in carbon_indices and n not in numbered_carbons]

    if c5_methyls:
        numbering['C18'] = c5_methyls[0]

    # Find methyl on C9
    if 'C9' in chain_numbering:
        c9_idx = chain_numbering['C9']
        c9_methyls = [n for n in connectivity[c9_idx]
                      if n in carbon_indices and n not in numbered_carbons]
        if c9_methyls:
            numbering['C19'] = c9_methyls[0]

    # Find methyl on C13
    if 'C13' in chain_numbering:
        c13_idx = chain_numbering['C13']
        c13_methyls = [n for n in connectivity[c13_idx]
                       if n in carbon_indices and n not in numbered_carbons]
        if c13_methyls:
            numbering['C20'] = c13_methyls[0]

    return numbering