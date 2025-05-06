#!/usr/bin/env python3
"""
This script reads opsin property data from CSV file and initializes PROTOS data structure.
It sets up the necessary datasets for the opsin analysis project.
"""

import os
import pandas as pd
import numpy as np
import json
import shutil
import pandas as pd
from pathlib import Path
import sys

project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root))

# Import PROTOS modules
from protos.io.paths.path_config import ProtosPaths, DataSource
from protos.processing.structure.struct_base_processor import CifBaseProcessor


def setup_registry_and_datasets(cp, paths, user_data_root, project_root):
    """Helper function to setup registry and load datasets"""
    # Get registry paths for both user and reference data
    user_registry_path = paths.get_registry_path("structure", DataSource.USER)
    ref_registry_path = paths.get_registry_path("structure", DataSource.REFERENCE)

    print(f"User registry path: {user_registry_path}")
    print(f"Reference registry path: {ref_registry_path}")

    # Convert string paths to Path objects if needed
    if isinstance(user_registry_path, str):
        user_registry_path = Path(user_registry_path)
    if isinstance(ref_registry_path, str):
        ref_registry_path = Path(ref_registry_path)

    # Create registry.json in both user and reference data if they don't exist
    for reg_path in [ref_registry_path, user_registry_path]:
        if not reg_path.exists():
            try:
                reg_path.parent.mkdir(parents=True, exist_ok=True)
                with open(reg_path, 'w') as f:
                    json.dump({"datasets": []}, f)
                print(f"Created empty registry at {reg_path}")
            except Exception as e:
                print(f"Warning: Could not create registry at {reg_path}: {str(e)}")

    # Import existing datasets
    # Check the legacy dataset location
    legacy_dataset_path = project_root / "protos" / "data" / "structure" / "structure_dataset" / "datasets.json"
    if legacy_dataset_path.exists():
        print(f"Found legacy dataset at: {legacy_dataset_path}")

        # Create structure directory in reference data using actual paths.ref_data_root
        ref_data_dir = Path(paths.ref_data_root) if paths.ref_data_root else Path(
            project_root) / "protos" / "src" / "protos" / "reference_data"
        struct_dataset_dir = ref_data_dir / "structure" / "structure_dataset"
        try:
            struct_dataset_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print(f"Warning: Could not create directory {struct_dataset_dir}: {str(e)}")

        # Copy the legacy dataset to the reference data directory
        target_dataset_path = struct_dataset_dir / "datasets.json"

        # Also create a copy in the user data directory
        user_struct_dataset_dir = Path(user_data_root) / "structure" / "structure_dataset"
        try:
            user_struct_dataset_dir.mkdir(parents=True, exist_ok=True)
            user_target_path = user_struct_dataset_dir / "datasets.json"
            shutil.copy2(legacy_dataset_path, user_target_path)
            print(f"Also copied dataset to user data: {user_target_path}")
        except Exception as e:
            print(f"Warning: Could not copy to user data: {str(e)}")
        try:
            shutil.copy2(legacy_dataset_path, target_dataset_path)
            print(f"Copied legacy dataset to: {target_dataset_path}")
        except Exception as e:
            print(f"Could not copy dataset file: {str(e)}")
            # If we can't copy, we'll just read from the original

        # Load datasets from the file (falling back to original if needed)
        try:
            dataset_file_to_read = target_dataset_path if target_dataset_path.exists() else legacy_dataset_path
            with open(dataset_file_to_read, 'r') as f:
                datasets_dict = json.load(f)
                print(f"Successfully loaded datasets from {dataset_file_to_read}")
        except Exception as e:
            print(f"Error reading datasets file: {str(e)}")
            datasets_dict = {}

        # Update the registry with the available datasets
        # First, create a standardized datasets file that matches the expected format
        print("Converting legacy dataset format to standardized format...")

        # Create a directory for standardized datasets
        std_dataset_dir = struct_dataset_dir / "standard"
        std_dataset_dir.mkdir(parents=True, exist_ok=True)

        # Create a dataset file for each dataset in the legacy format
        registry_data = {"datasets": []}
        for dataset_name, pdb_ids in datasets_dict.items():
            # Create a standardized dataset file
            std_dataset_path = std_dataset_dir / f"{dataset_name}.json"

            # Format the dataset as a dictionary with metadata and pdb_ids array
            std_dataset = {
                "id": dataset_name,
                "name": dataset_name,
                "description": f"Dataset containing {len(pdb_ids)} structures",
                "type": "structure",
                "pdb_ids": pdb_ids
            }

            # Write the standardized dataset file
            with open(std_dataset_path, 'w') as f:
                json.dump(std_dataset, f, indent=2)
            print(f"Created standardized dataset file: {std_dataset_path}")

            # Add to registry
            registry_data["datasets"].append({
                "id": dataset_name,
                "path": str(std_dataset_path),
                "type": "json"
            })

        # Save the updated registry to both reference and user data
        for reg_path in [ref_registry_path, user_registry_path]:
            try:
                with open(reg_path, 'w') as f:
                    json.dump(registry_data, f, indent=2)
                print(f"Updated registry at {reg_path}")
            except Exception as e:
                print(f"Warning: Could not update registry at {reg_path}: {str(e)}")

        # Also directly update the processor's registry
        if hasattr(cp, "registry"):
            cp.registry = registry_data
            print("Updated processor registry directly")

        # Reinitialize the processor to pick up the new registry
        try:
            cp.__init__(
                name="opsin_processor",
                data_root=str(user_data_root),  # Convert Path to string for processor
                processor_data_dir="structure",  # Standard subdirectory
                preload=False  # Don't load PDB IDs on init
            )
            print("Reinitialized processor to pick up new registry")
        except Exception as e:
            print(f"Error reinitializing processor: {str(e)}")


def setup_processor(user_data_root):
    """
    Set up the CifBaseProcessor with proper path resolution
    """
    print("=== Setting up CifBaseProcessor with Modern Path Resolution ===")

    # Ensure user_data_root is a Path object and is absolute
    if isinstance(user_data_root, str):
        user_data_root = Path(user_data_root)
    user_data_root = user_data_root.resolve()
    
    print(f"Using user_data_root: {user_data_root}")
    
    # Let ProtosPaths find the reference data root by passing None
    # This will use the appropriate default or environment variable
    ref_data_root = None
    
    # Create an instance of ProtosPaths with proper settings
    paths = ProtosPaths(
        user_data_root=str(user_data_root),  # Convert to string as expected by ProtosPaths
        ref_data_root=ref_data_root,
        create_dirs=True,  # Create directories that don't exist
        validate=True      # Validate directory structure (shows warnings for missing directories)
    )

    print(f"ProtosPaths initialized:")
    print(f"  User data root: {paths.user_data_root}")
    print(f"  Reference data root: {paths.ref_data_root}")

    # Initialize the processor with modern path parameters
    # Note: BaseProcessor internally creates its own ProtosPaths instance
    # using the user_data_root we provide
    cp = CifBaseProcessor(
        name="opsin_processor",
        data_root=str(user_data_root),  # Convert to string as expected by CifBaseProcessor
        processor_data_dir="structure",  # Standard subdirectory
        preload=False  # Don't load PDB IDs on init
    )

    # Set up registry and datasets
    setup_registry_and_datasets(cp, paths, user_data_root, project_root)

    return cp, paths

def read_opsin_csv(csv_path=None):
    """Read the opsin property CSV file."""
    # Default path to the CSV file if not provided
    if csv_path is None:
        csv_path = Path(project_root) / "property" / "mo_exp.csv"
    
    # Ensure csv_path is a Path object
    if not isinstance(csv_path, Path):
        csv_path = Path(csv_path)
    
    # Check if file exists
    if not csv_path.exists():
        # Try alternate locations
        alternate_paths = [
            project_root / "projects" / "opsin_analysis" / "property" / "mo_exp.csv",
            project_root / "mo_exp.csv"
        ]
        
        for alt_path in alternate_paths:
            if alt_path.exists():
                csv_path = alt_path
                print(f"Found CSV file at alternate location: {csv_path}")
                break
        else:
            print(f"Warning: CSV file not found at {csv_path} or any alternate locations")
            return pd.DataFrame()
    
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Print columns and some basic info
    print(f"CSV file columns: {df.columns.tolist()}")
    print(f"Total entries: {len(df)}")
    
    return df


def process_local_dataset(cp, dataset_name, dataset_info, user_data_root):
    """Process structures from local files"""

    # Convert relative path to absolute
    source_path = Path(dataset_info['source_path'])
    if not source_path.is_absolute():
        source_path = Path(project_root) / source_path
    
    if not source_path.exists():
        print(f"Source path not found: {source_path}")
        return

    # Get all CIF files
    cif_files = list(source_path.glob("*.cif"))
    if not cif_files:
        print(f"No CIF files found in {source_path}")
        return

    print(f"Found {len(cif_files)} CIF files for {dataset_name}")

    # Create the destination directory in the PROTOS structure
    mmcif_dir = user_data_root / "structure" / "mmcif"
    mmcif_dir.mkdir(exist_ok=True, parents=True)

    # Copy files to the PROTOS structure directory
    pdb_ids = []
    successful_copies = 0
    
    for cif_file in cif_files:
        # Generate a PROTOS-compatible ID - prevent scientific notation conversion
        pdb_id = cif_file.stem
        
        # Check if PDB ID could be interpreted as scientific notation
        if any(c.lower() == 'e' for c in pdb_id):
            print(f"Warning: PDB ID '{pdb_id}' contains 'e' and may be interpreted as scientific notation")
        
        # Copy the file to the PROTOS directory - ensure exact filename is preserved
        dest_file = mmcif_dir / cif_file.name
        try:
            # Only copy if target doesn't exist or is older
            if not dest_file.exists() or (cif_file.stat().st_mtime > dest_file.stat().st_mtime):
                shutil.copy2(cif_file, dest_file)
                print(f"Copied {cif_file.name} to {dest_file}")
                successful_copies += 1
            else:
                print(f"Skipped {cif_file.name} (already exists and is up to date)")
            
            # Only add to the list if copy was successful or file already exists
            pdb_ids.append(pdb_id)
        except Exception as e:
            print(f"Error copying {cif_file.name}: {e}")

    print(f"Copied {successful_copies} files, skipped {len(cif_files) - successful_copies} files")
    
    if not pdb_ids:
        print(f"No structures were successfully copied for {dataset_name}")
        return
        
    # Load structures into the processor
    try:
        # Reset processor to ensure a clean state
        cp.reset_data()
        
        # Ensure PDB IDs are preserved as strings (not converted to scientific notation)
        pdb_ids_str = [str(pdb_id) for pdb_id in pdb_ids]
        
        # Check for potential scientific notation issues
        for pdb_id in pdb_ids_str:
            if any(c.lower() == 'e' for c in pdb_id):
                print(f"Warning: PDB ID '{pdb_id}' may be interpreted as scientific notation")
        
        # Load the structures
        cp.load_structures(pdb_ids_str)
        
        # Check if any structures were loaded
        if len(cp.pdb_ids) > 0:
            # Create dataset with standardized metadata
            metadata = {
                "description": dataset_info.get('description', ''),
                "source": dataset_info.get('source_path', ''),
                "creation_date": pd.Timestamp.now().strftime('%Y-%m-%d'),
                "structures_count": len(cp.pdb_ids)
            }
            
            # Create the dataset with metadata
            cp.create_dataset(
                dataset_id=dataset_name,
                name=dataset_name,
                description=metadata['description'],
                content=cp.pdb_ids,
                metadata=metadata
            )
            
            # Create a standardized dataset file as well
            dataset_dir = user_data_root / "structure" / "structure_dataset" / "standard"
            dataset_dir.mkdir(exist_ok=True, parents=True)
            
            std_dataset_path = dataset_dir / f"{dataset_name}.json"
            std_dataset = {
                "id": dataset_name,
                "name": dataset_name,
                "description": metadata['description'],
                "type": "structure",
                "pdb_ids": pdb_ids,
                "metadata": metadata
            }
            
            with open(std_dataset_path, 'w') as f:
                json.dump(std_dataset, f, indent=2)
                
            print(f"Successfully saved dataset {dataset_name} with {len(cp.pdb_ids)} structures")
            print(f"Created standardized dataset file: {std_dataset_path}")
        else:
            print(f"No structures were successfully loaded for {dataset_name}")
    except Exception as e:
        print(f"Error processing dataset {dataset_name}: {e}")


def process_downloadable_dataset(cp, dataset_name, dataset_info):
    """Process structures that need to be downloaded from PDB"""

    pdb_ids = dataset_info.get('pdb_ids', [])
    if not pdb_ids:
        print(f"No PDB IDs provided for {dataset_name}")
        return

    print(f"Downloading {len(pdb_ids)} structures for {dataset_name}")

    # Reset processor to ensure clean state
    cp.reset_data()

    # Track successful downloads
    successful_ids = []
    
    # Download structures in batches to avoid timeout issues
    batch_size = 10
    total_batches = (len(pdb_ids) + batch_size - 1) // batch_size
    
    for i in range(0, len(pdb_ids), batch_size):
        batch = pdb_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        print(f"Processing batch {batch_num}/{total_batches} [{len(batch)} structures]")

        try:
            # Create a temporary dataset with these IDs
            temp_dataset_name = f"{dataset_name}_temp_{batch_num}"
            
            # Create dataset with batch IDs
            cp.reset_data()
            for pdb_id in batch:
                if pdb_id not in cp.pdb_ids:
                    cp.pdb_ids.append(pdb_id)
            
            # Create temporary dataset
            cp.create_dataset(
                dataset_id=temp_dataset_name,
                name=temp_dataset_name,
                description="Temporary batch for downloading",
                content=batch
            )
            
            # Download missing CIFs for this dataset 
            downloaded = cp.check_and_download_missing_cifs(temp_dataset_name)
            
            # Reload the batch to make sure we have the actual data
            cp.reset_data()
            
            # Ensure batch PDB IDs are preserved as strings (not converted to scientific notation)
            batch_str = [str(pdb_id) for pdb_id in batch]
            
            # Check for potential scientific notation issues
            for pdb_id in batch_str:
                if any(c.lower() == 'e' for c in pdb_id):
                    print(f"Warning: PDB ID '{pdb_id}' may be interpreted as scientific notation")
            
            cp.load_structures(batch_str)
            
            # Check which ones were successfully loaded
            batch_successful = [pdb_id for pdb_id in batch if pdb_id in cp.pdb_ids]
            successful_ids.extend(batch_successful)
            
            if batch_successful:
                print(f"Successfully downloaded {len(batch_successful)}/{len(batch)} structures in batch {batch_num}")
            else:
                print(f"Failed to download structures in batch {batch_num}")
                
        except Exception as e:
            print(f"Error downloading batch {batch_num}: {e}")

    # After downloading all structures, save the dataset if any were successful
    if successful_ids:
        try:
            # Create dataset with standardized metadata
            metadata = {
                "description": dataset_info.get('description', ''),
                "source": "RCSB PDB",
                "creation_date": pd.Timestamp.now().strftime('%Y-%m-%d'),
                "structures_count": len(successful_ids),
                "original_count": len(pdb_ids),
                "success_rate": f"{len(successful_ids)}/{len(pdb_ids)} ({100 * len(successful_ids) / len(pdb_ids):.1f}%)"
            }
            
            # Create the dataset with metadata
            cp.create_dataset(
                dataset_id=dataset_name,
                name=dataset_name,
                description=metadata['description'],
                content=successful_ids,
                metadata=metadata
            )
            
            # Create a standardized dataset file as well
            user_data_root = Path(cp.data_root)
            dataset_dir = user_data_root / "structure" / "structure_dataset" / "standard"
            dataset_dir.mkdir(exist_ok=True, parents=True)
            
            std_dataset_path = dataset_dir / f"{dataset_name}.json"
            std_dataset = {
                "id": dataset_name,
                "name": dataset_name,
                "description": metadata['description'],
                "type": "structure",
                "pdb_ids": successful_ids,
                "metadata": metadata
            }
            
            with open(std_dataset_path, 'w') as f:
                json.dump(std_dataset, f, indent=2)
                
            print(f"Successfully saved dataset {dataset_name} with {len(successful_ids)} structures")
            print(f"Created standardized dataset file: {std_dataset_path}")
        except Exception as e:
            print(f"Error saving dataset {dataset_name}: {e}")
    else:
        print(f"No structures were successfully downloaded for {dataset_name}")
        print(f"Dataset {dataset_name} not created")


def validate_datasets(cp, datasets, user_data_root):
    """Validate that all datasets were created properly"""

    print("\n=== Validating Datasets ===")
    
    # List available datasets
    try:
        # Get the list of dataset names from registry
        if hasattr(cp, "registry") and "datasets" in cp.registry:
            available_datasets = [d["id"] for d in cp.registry["datasets"] if "id" in d]
        else:
            available_datasets = []
        
        print(f"Available datasets: {available_datasets}")
    except Exception as e:
        print(f"Error listing datasets: {e}")
        return

    for dataset_name in datasets:
        if dataset_name in available_datasets:
            # Try loading the dataset to verify it's accessible
            cp.reset_data()
            try:
                # Load the dataset
                cp.load_dataset(dataset_name)
                print(f"Dataset {dataset_name} loaded successfully with {len(cp.pdb_ids)} structures")
            except Exception as e:
                print(f"Error loading dataset {dataset_name}: {e}")
        else:
            print(f"Dataset {dataset_name} not found in available datasets")
            
            # Check if the standard dataset file exists directly
            std_dataset_path = user_data_root / "structure" / "structure_dataset" / "standard" / f"{dataset_name}.json"
            if std_dataset_path.exists():
                print(f"Standard dataset file exists at {std_dataset_path} but it's not registered in the processor")

def main():
    """
    Main function for setting up opsin analysis datasets
    """
    # Initialize paths and data - use absolute paths
    project_root = Path(__file__).resolve().parent
    user_data_root = project_root / "data"
    
    # Ensure user_data_root exists
    user_data_root.mkdir(parents=True, exist_ok=True)
    
    # Look for mo_exp.csv in multiple possible locations
    csv_path = project_root / "property" / "mo_exp.csv"
    if not csv_path.exists():
        alternate_paths = [
            project_root / "projects" / "opsin_analysis" / "property" / "mo_exp.csv",
            project_root / "mo_exp.csv"
        ]
        
        for alt_path in alternate_paths:
            if alt_path.exists():
                csv_path = alt_path
                print(f"Found CSV file at: {csv_path}")
                break
        else:
            print(f"Warning: CSV file not found at expected locations. Will try with provided path.")
            csv_path = Path("projects/opsin_analysis/property/mo_exp.csv")

    # Setup processor with proper error handling
    try:
        cp, paths = setup_processor(user_data_root)
        print("\n=== CifBaseProcessor Setup Complete ===")
        print(f"Processor type: {type(cp).__name__}")
        print(f"User data root: {paths.user_data_root}")
        print(f"Reference data root: {paths.ref_data_root}")
    except Exception as e:
        print(f"Failed to initialize processor: {e}")
        return

    # Read property data using the improved function
    try:
        df = read_opsin_csv(csv_path)
        print(f"\n=== Property Data ===")
        print(f"Unique PDB IDs: {len(df['pdb_id'].unique())}")
    except Exception as e:
        print(f"Error reading property data: {e}")
        return

    # Ensure required directories exist
    mmcif_dir = user_data_root / "structure" / "mmcif"
    dataset_dir = user_data_root / "structure" / "structure_dataset"
    standard_dir = dataset_dir / "standard"
    
    for directory in [mmcif_dir, dataset_dir, standard_dir]:
        directory.mkdir(exist_ok=True, parents=True)
        print(f"Ensured directory exists: {directory}")

    # Process all four datasets with proper error handling
    datasets = {
        'mo_exp': {
            'description': 'Experimental microbial opsin structures',
            'pdb_ids': df['pdb_id'].dropna().unique().tolist(),
            'source_type': 'download'
        },
        'mo_pred': {
            'description': 'Predicted microbial opsin structures',
            'source_path': 'structures/mo_pred',  # Relative to project_root
            'source_type': 'local'
        },
        'hideaki_exp': {
            'description': 'Experimental structures from Hideaki dataset',
            'source_path': 'structures/hideaki_exp',  # Relative to project_root
            'source_type': 'local'
        },
        'hideaki_pred': {
            'description': 'Predicted structures from Hideaki dataset',
            'source_path': 'structures/hideaki_pred',  # Relative to project_root
            'source_type': 'local'
        }
    }

    # Process each dataset
    for name, info in datasets.items():
        print(f"\n=== Processing {name} dataset ===")

        # Clear processor data for new dataset
        cp.reset_data()

        # Handle based on source type
        if info['source_type'] == 'download':
            process_downloadable_dataset(cp, name, info)
        else:
            process_local_dataset(cp, name, info, user_data_root)

    # Final validation
    validate_datasets(cp, datasets.keys(), user_data_root)

    print("\nAll datasets have been successfully prepared!")
if __name__ == "__main__":
    main()