from src.data_processing import load_opsin_property_data
import pickle
import os
import sys
from pathlib import Path
import pandas as pd
import traceback


current_file_path = Path(__file__).resolve()
PROJECT_DIR = current_file_path.parent

CACHE_DIR = PROJECT_DIR / 'opsin_output' / 'cache'
GRN_TABLES_DIR = PROJECT_DIR / 'opsin_output' / 'opsin_grn_tables'

project_dir = Path(__file__).resolve().parent
src_dir = project_dir / 'src'
if str(src_dir) not in sys.path:
    sys.path.append(str(src_dir))


def load_cached_data(cache_path, description="data"):
    if os.path.exists(cache_path):
        print(f"[INFO] Loading {description} from cache: {cache_path}")
        try:
            with open(cache_path, 'rb') as f:
                result = pickle.load(f)
            print(f"[INFO] Successfully loaded {description}")
            return result
        except Exception as e:
            print(f"[ERROR] Error loading {description} from cache '{cache_path}': {e}")
            traceback.print_exc()
    else:
        print(f"[WARN] Cache file not found: {cache_path}")
    return None


# Check property data and structure mapping
print("\n=== Property Data and Structure Mapping ===")

# Load processed structures to pass to property loading
processed_structures = {}
processed_path = 'opsin_output/cache/processed_structures_A.pkl'
if os.path.exists(processed_path):
    with open(processed_path, 'rb') as f:
        data = pickle.load(f)
        if 'processed_structures' in data:
            processed_structures = data['processed_structures']

print(len(processed_structures))
print(processed_structures.keys())
print(list(processed_structures['4PXK'].keys())[:10])
print(processed_structures['4PXK']['structure_type'])
print(processed_structures['4PXK']['base_name'])
print(processed_structures['HsHR_model_0']['base_name'])
print(processed_structures['R2ACR_J315_refine8']['base_name'])

property_csv_path = project_dir / 'property' / 'mo_exp.csv'
if property_csv_path.exists():
    property_data = load_opsin_property_data(property_csv_path, processed_structures)

print(property_data['properties']['3UG9'])


def load_grn_tables_data():  # Matching notebook name
    grn_tables_pkl = GRN_TABLES_DIR / 'grn_tables.pkl'  # Corrected path
    return load_cached_data(grn_tables_pkl, "GRN tables data")


grn_tables = load_grn_tables_data()
print(grn_tables.keys())
print(grn_tables['residue_table'].head(2))
print(grn_tables['helix_pivot_columns'])
