
from pathlib import Path
import os
from protos.io.paths.path_config import ProtosPaths, DataSource
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Get absolute path to user data directory
user_data_root = Path(__file__).resolve().parent / "data"
abs_user_data_root = str(user_data_root.absolute())

# Set environment variable to force ProtosPaths to use this directory
os.environ["PROTOS_DATA_ROOT"] = abs_user_data_root
os.environ["PROTOS_REF_DATA_ROOT"] = abs_user_data_root  # Optional, to force reference to be same as user data

# Create paths explicitly forcing USER data source
paths = ProtosPaths(
    user_data_root=abs_user_data_root,
    ref_data_root=abs_user_data_root,  # Force reference data to be the same as user data
    create_dirs=True
)

print(f"User data root: {paths.user_data_root}")
print(f"Reference data root: {paths.ref_data_root}")

# Initialize the processor with explicit user data path
cp = CifBaseProcessor(
    name="cif_processor",
    data_root=abs_user_data_root,
    processor_data_dir="structure",
)

# Override key CifBaseProcessor paths to ensure consistency
cp.path_structure_dir = os.path.join(abs_user_data_root, "structure", "mmcif")
cp.path_dataset_dir = os.path.join(abs_user_data_root, "structure", "structure_dataset")
cp.path_alignment_dir = os.path.join(abs_user_data_root, "structure", "alignments")

# Check standard directories exist
os.makedirs(cp.path_structure_dir, exist_ok=True)
os.makedirs(cp.path_dataset_dir, exist_ok=True)
os.makedirs(cp.path_alignment_dir, exist_ok=True)

# List datasets
print(f"Processor data_root: {cp.data_root}")
print(f"Processor dataset dir: {cp.path_dataset_dir}")
datasets = cp.list_datasets()
print("Datasets:", datasets)

# Check registry exists in user data path
registry_path = Path(abs_user_data_root) / "structure" / "registry.json"
print(f"Registry path exists: {registry_path.exists()}")

# If needed, manually scan for datasets in structure_dataset directory
structure_dataset_dir = Path(abs_user_data_root) / "structure" / "structure_dataset" / "standard"
if structure_dataset_dir.exists():
    print(f"Available dataset files:")
    for file in structure_dataset_dir.glob("*.json"):
        print(f"  {file.name}")

cp.load_dataset('hideaki_exp')
cp.load_dataset('hideaki_pred')
cp.load_dataset('mo_exp')
cp.load_dataset('mo_pred')
