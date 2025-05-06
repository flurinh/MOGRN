# CLAUDE.md - Protos Framework Reference

## ProtosPaths System

ProtosPaths handles path management throughout the Protos framework:

```python
# Key classes in protos.io.paths
from protos.io.paths.path_config import ProtosPaths, DataSource

# Initialization
paths = ProtosPaths(
    user_data_root=None,  # Uses env var or defaults to 'data'
    ref_data_root=None,   # Uses package resources
    create_dirs=True,     # Creates directories if missing
    validate=True         # Validates directory structure
)

# Core functions
structure_path = paths.get_processor_path("structure")
abs_path = paths.resolve_path("dataset.json", source=DataSource.USER)
exists, source = paths.exists("structure/1abc.cif")

# Global helpers (from protos.io.paths)
structure_file = get_structure_path("1abc")
dataset_file = get_dataset_path("my_dataset", processor_type="structure")
```

## DatasetManager

DatasetManager provides standardized dataset operations:

```python
# In protos.core.dataset_manager
manager = DatasetManager(
    processor_type="structure",  # Processor type
    paths=paths                  # Path resolver instance
)

# Dataset operations
dataset = manager.create_dataset(
    dataset_id="my_dataset",
    name="My Structure Dataset",
    description="A collection of structure IDs",
    content=["1abc", "2xyz", "3def"],
    metadata={"source": "PDB"}
)

dataset = manager.load_dataset("my_dataset")
structures = dataset.content
manager.save_dataset(dataset)
datasets = manager.list_datasets()
manager.delete_dataset("my_dataset")
```

## BaseProcessor

BaseProcessor is the foundation for all processor types:

```python
# In protos.core.base_processor
class StructureProcessor(BaseProcessor):
    def __init__(self, name, data_root=None, config=None):
        super().__init__(
            name=name,
            data_root=data_root,
            processor_data_dir="structure", # Directory name
            config=config
        )

# Initialization pattern for processors
processor = StructureProcessor(
    name="my_processor",     # Identifier
    data_root="~/data",      # Custom data location (optional)
    config={"param": "val"}  # Configuration values
)

# Common methods
data = processor.load_data("dataset_id", file_format="csv")
processor.save_data("results", data=processed_data, file_format="json")
dataset = processor.create_standard_dataset(
    dataset_id="std_dataset",
    name="Standard Dataset",
    content=["item1", "item2"]
)
```

## Integration Flow

1. BaseProcessor creates ProtosPaths instance during initialization
2. BaseProcessor initializes DatasetManager with paths and processor type
3. When loading/saving:
   - ProtosPaths resolves file paths
   - DatasetManager handles dataset operations
   - BaseProcessor manages format conversion and error handling