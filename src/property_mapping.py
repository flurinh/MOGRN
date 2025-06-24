"""
Unified property mapping module for consistent structure-to-property association.
Handles 4 distinct dataset types and their relationships:
1. Standard experimental structures (PDB IDs)
2. Standard predicted structures (_model_0)
3. Hideaki experimental structures (_refine)
4. Hideaki predicted structures (_refine_model_0)
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


class PropertyMapper:
    """Handles mapping between structure IDs and property data."""
    
    def __init__(self, property_csv_path: Path):
        """Initialize with property CSV file."""
        self.property_df = pd.read_csv(property_csv_path)
        self.property_cache = {}
        self.exp_pred_pairs = {}  # Maps experimental -> predicted
        self._build_mapping()
    
    def identify_structure_type(self, structure_id: str) -> Tuple[str, str]:
        """
        Identify structure type and extract base name.
        
        Returns:
            (structure_type, base_name)
            
        Structure types:
            - 'standard_exp': Standard experimental (PDB ID)
            - 'standard_pred': Standard predicted (_model_0)
            - 'hideaki_exp': Hideaki experimental (_refine)
            - 'hideaki_pred': Hideaki predicted (_refine_model_0)
        """
        # Check Hideaki predicted first (most specific pattern)
        hideaki_pred_pattern = r'^(.+)_J\d+_refine\d+_model_0$'
        match = re.match(hideaki_pred_pattern, structure_id)
        if match:
            return 'hideaki_pred', match.group(1)
        
        # Check Hideaki experimental
        hideaki_exp_pattern = r'^(.+)_J\d+_refine\d+$'
        match = re.match(hideaki_exp_pattern, structure_id)
        if match:
            return 'hideaki_exp', match.group(1)
        
        # Check standard predicted (including _smile variants)
        if structure_id.endswith('_model_0') or structure_id.endswith('_smile_model_0'):
            base = structure_id
            for suffix in ['_smile_model_0', '_model_0']:
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
            return 'standard_pred', base
        
        # Default to standard experimental
        return 'standard_exp', structure_id
    
    def clean_name(self, name: str) -> str:
        """Clean structure name for consistent matching."""
        if not name:
            return name
        cleaned = str(name).strip()
        cleaned = cleaned.replace('.', '')
        cleaned = cleaned.replace('+', '')
        cleaned = cleaned.replace('-', '_')
        return cleaned
    
    def _build_mapping(self):
        """Build comprehensive mapping for all 4 dataset types."""
        for idx, row in self.property_df.iterrows():
            # Extract fields
            pdb_id = str(row['pdb_id']).strip() if pd.notna(row['pdb_id']) else None
            short_name = str(row['short_name']).strip() if pd.notna(row['short_name']) else None
            opsin_name = str(row['opsin name']).strip() if pd.notna(row['opsin name']) else None
            
            # Build property data
            property_data = {
                'domain': str(row['Rhodopsin Type (Microbial)']).strip() 
                          if pd.notna(row['Rhodopsin Type (Microbial)']) else 'Unknown',
                'molecular_function': str(row['molecular_function']).strip().replace('?', '')
                                    if pd.notna(row['molecular_function']) else 'Unknown',
                'experimentally_determined': bool(row.get('experimentally_determined', False))
                                           if 'experimentally_determined' in row else False,
                'short_name': short_name,
                'opsin_name': opsin_name,
                'pdb_id': pdb_id
            }
            
            # Handle standard structures (experimental + predicted)
            if pdb_id:
                # Map experimental structure
                self.property_cache[pdb_id] = property_data.copy()
                self.property_cache[pdb_id.upper()] = property_data.copy()
                self.property_cache[pdb_id.lower()] = property_data.copy()
                
                # Map corresponding predicted structure
                if short_name:
                    cleaned_short = self.clean_name(short_name)
                    pred_id = f"{cleaned_short}_model_0"
                    self.property_cache[pred_id] = property_data.copy()
                    
                    # Create experimental-predicted pair
                    self.exp_pred_pairs[pdb_id] = pred_id
                    
                    # Also handle _smile variant for compatibility
                    smile_pred_id = f"{cleaned_short}_smile_model_0"
                    self.property_cache[smile_pred_id] = property_data.copy()
            
            # Handle Hideaki structures (use short_name as base)
            elif short_name:
                cleaned_name = self.clean_name(short_name)
                
                # This entry represents both Hideaki experimental and predicted
                # The actual structure IDs will have _J###_refine# patterns
                # Store under the base name for lookup
                self.property_cache[short_name] = property_data.copy()
                self.property_cache[cleaned_name] = property_data.copy()
                
                # Also store standard predicted variant
                self.property_cache[f"{cleaned_name}_model_0"] = property_data.copy()
    
    def get_properties(self, structure_id: str) -> Optional[Dict]:
        """Get properties for any structure ID."""
        # Direct lookup
        if structure_id in self.property_cache:
            return self.property_cache[structure_id]
        
        # Identify structure type and base name
        struct_type, base_name = self.identify_structure_type(structure_id)
        
        # Try base name
        if base_name in self.property_cache:
            return self.property_cache[base_name]
        
        # Try cleaned base name
        cleaned_base = self.clean_name(base_name)
        if cleaned_base in self.property_cache:
            return self.property_cache[cleaned_base]
        
        # Special handling for case sensitivity
        if base_name.upper() in self.property_cache:
            return self.property_cache[base_name.upper()]
        
        # Try without _smile
        no_smile = structure_id.replace('_smile', '')
        if no_smile in self.property_cache:
            return self.property_cache[no_smile]
        
        return None
    
    def get_structure_pairs(self) -> Dict[str, List[Tuple[str, str]]]:
        """
        Get all experimental-predicted structure pairs organized by type.
        
        Returns dict with keys:
            - 'standard': List of (exp_id, pred_id) tuples
            - 'hideaki': List of (exp_id, pred_id) tuples
        """
        return {
            'standard': list(self.exp_pred_pairs.items()),
            'hideaki': []  # Will be populated when we scan actual structures
        }
    
    def find_hideaki_pairs(self, all_structure_ids: List[str]) -> List[Tuple[str, str]]:
        """
        Find Hideaki experimental-predicted pairs from structure IDs.
        
        Returns list of (hideaki_exp, hideaki_pred) tuples.
        """
        hideaki_pairs = []
        hideaki_exp_ids = []
        
        # First find all Hideaki experimental structures
        for struct_id in all_structure_ids:
            struct_type, base_name = self.identify_structure_type(struct_id)
            if struct_type == 'hideaki_exp':
                hideaki_exp_ids.append(struct_id)
        
        # For each experimental, find its predicted counterpart
        for exp_id in hideaki_exp_ids:
            pred_id = f"{exp_id}_model_0"
            if pred_id in all_structure_ids:
                hideaki_pairs.append((exp_id, pred_id))
        
        return hideaki_pairs
    
    def update_processed_structures(self, processed_structures: Dict) -> Dict:
        """Update processed structures with property data."""
        updated_count = 0
        missing_properties = []
        
        # Also find Hideaki pairs
        all_ids = list(processed_structures.keys())
        hideaki_pairs = self.find_hideaki_pairs(all_ids)
        
        for struct_id, struct_data in processed_structures.items():
            properties = self.get_properties(struct_id)
            
            if properties:
                if 'properties' not in struct_data:
                    struct_data['properties'] = {}
                struct_data['properties'].update(properties)
                updated_count += 1
            else:
                missing_properties.append(struct_id)
        
        print(f"Updated properties for {updated_count} structures")
        if missing_properties:
            print(f"Missing properties for {len(missing_properties)} structures:")
            print(f"Examples: {missing_properties[:5]}")
        
        if hideaki_pairs:
            print(f"\nFound {len(hideaki_pairs)} Hideaki experimental-predicted pairs")
        
        return processed_structures
    
    def get_all_experimental_predicted_pairs(self, all_structure_ids: List[str]) -> Dict[str, str]:
        """
        Get complete mapping of experimental to predicted structures.
        Includes both standard and Hideaki pairs.
        """
        all_pairs = self.exp_pred_pairs.copy()
        
        # Add Hideaki pairs
        hideaki_pairs = self.find_hideaki_pairs(all_structure_ids)
        for exp_id, pred_id in hideaki_pairs:
            all_pairs[exp_id] = pred_id
        
        return all_pairs
    
    def get_all_properties(self) -> Dict[str, Dict]:
        """Get all property mappings."""
        return self.property_cache.copy()
    
    def validate_mapping(self, processed_structures: Dict) -> Dict:
        """Validate property mapping and structure relationships."""
        all_ids = list(processed_structures.keys())
        stats = {
            'total_structures': len(processed_structures),
            'structures_with_properties': 0,
            'structures_without_properties': 0,
            'standard_exp': 0,
            'standard_pred': 0,
            'hideaki_exp': 0,
            'hideaki_pred': 0,
            'missing_properties': [],
            'structure_types': {}
        }
        
        for struct_id, struct_data in processed_structures.items():
            # Check structure type
            struct_type, base_name = self.identify_structure_type(struct_id)
            stats['structure_types'][struct_id] = struct_type
            stats[struct_type] = stats.get(struct_type, 0) + 1
            
            # Check properties
            has_properties = False
            if 'properties' in struct_data and struct_data['properties']:
                has_properties = True
            elif self.get_properties(struct_id):
                has_properties = True
            
            if has_properties:
                stats['structures_with_properties'] += 1
            else:
                stats['structures_without_properties'] += 1
                stats['missing_properties'].append(struct_id)
        
        # Find all pairs
        hideaki_pairs = self.find_hideaki_pairs(all_ids)
        stats['standard_pairs'] = len(self.exp_pred_pairs)
        stats['hideaki_pairs'] = len(hideaki_pairs)
        stats['total_pairs'] = stats['standard_pairs'] + stats['hideaki_pairs']
        
        return stats


def create_unified_property_mapper(property_csv_path: Path) -> PropertyMapper:
    """Create and return a unified property mapper instance."""
    return PropertyMapper(property_csv_path)