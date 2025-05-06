# MOGRN Project Guide

## Step 1: Check Property Data and Prepare Protos Configuration

**Status:** ✅ COMPLETED

This step sets up the proper data structure and Protos framework configuration:

### Path Configuration System:
- **ProtosPaths System:** The project uses the ProtosPaths system from the Protos framework to manage all file paths
- **Data Organization:** 
  - `user_data_root`: Primary location for user-writable data (defaults to "./data")
  - `ref_data_root`: Reference data that is packaged with the framework
  - Processor-specific subdirectories (structure, grn, sequence, etc.)

### How Path Resolution Works:
- **Absolute Paths:** All working paths are resolved to absolute paths internally
- **Data Sources:** The system distinguishes between:
  - `DataSource.USER`: User-writable data
  - `DataSource.REFERENCE`: Read-only reference data
  - `DataSource.AUTO`: Automatically chooses the appropriate source
- **Directory Structure:** Standard directories are created for each processor type

### Using the System:
```python
# Core path resolution
paths = ProtosPaths(user_data_root="./data", create_dirs=True)
structure_path = paths.get_processor_path("structure")  # Get structure processor dir
registry_path = paths.get_registry_path("structure")    # Get registry file path

# Dataset management
cp = CifBaseProcessor(name="processor_name", data_root="./data")
cp.create_dataset(dataset_id="my_dataset", content=pdb_ids)
cp.load_dataset("my_dataset")  # Load structures from the dataset
```

The `prepare_data.py` script:
1. Initializes the ProtosPaths system with the proper root directories
2. Sets up the CifBaseProcessor with these paths
3. Processes opsin datasets and registers them in the system
4. Creates standardized dataset files in the proper locations

## Step 2: Create YAML Files

**Status:** 🔄 TODO

## Step 3: Set Up Protos Data Folder

**Status:** 🔄 TODO

## Step 4: Preprocess Structures/Dataset

**Status:** 🔄 TODO

## Step 5: Analysis Workflow

**Status:** 🔄 TODO

## Step 6: Visualizations

**Status:** 🔄 TODO

## Step 7: Assign GRNs to Dataset Using Protos

**Status:** 🔄 TODO