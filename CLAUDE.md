# CLAUDE.MD - QUICK REFERENCE

## 🔴 MANDATORY SESSION START PROCEDURE
**EVERY new session MUST begin with:**
1. **READ** `/claude/execution.txt` - Understand system setup and activation procedures
2. **READ** `/claude/content.txt` - Get current project structure
3. **READ** `/claude/logs.txt` - Check last 20-30 entries for recent context (use tail -n 30)
4. **READ** `/claude/todo.txt` - Review project-level task trees and pending work
5. **READ** `/claude/csv_tables.txt` - Understand CSV handling rules and context-aware operations

**This is NOT optional - Do this BEFORE any other action**

## 🚨 CRITICAL RULES - CHECK EVERY TIME
1. **LOG EVERY FILE CHANGE** → Append to `/claude/logs.txt` → `[YYYY-MM-DD HH:MM] [CREATE/MODIFY/DELETE] file - message`
2. **UPDATE PROJECT STRUCTURE** → `/claude/content.txt` → After any file creation/deletion
3. **TEMP FILES** → Prefix with `temp_` or `dev_`
4. **DEPRECATED FILES** → Prefix with `deprecated_`
5. **BEFORE ANY FILE OP** → Check `/claude/content.txt` first
6. **NO SUDO COMMANDS** → Claude cannot run sudo commands. User must manually run:
   - `sudo service postgresql start` (start database)
   - Any other system-level commands requiring sudo

## 🎯 PROJECT QUICK INFO

## 📁 KEY PATHS
```
/claude/logs.txt          → File change log (MANDATORY)
/claude/content.txt       → Current architecture
/claude/todo.txt          → Project task trees (XML format)
/claude/csv_tables.txt    → CSV handling rules and context
/claude/csv_helper.py     → CSV operations helper script
```

## 🔧 COMMON COMMANDS

### CSV Operations (Memory-Efficient)
```bash
# First understand table structure
python /claude/csv_helper.py schema output/rmsd_matrix.csv

# Get row by index as dictionary
python /claude/csv_helper.py get-row output/msa_table_grn.csv --index "PDB_001"

# Get specific cell value
python /claude/csv_helper.py get-cell output/distance_table_grn.csv --index "PDB_001" --column "3.50"

# Search indices
python /claude/csv_helper.py search-index output/protein_summary.csv --pattern ".*_exp$"
```

## 🔍 LOG SEARCH PATTERNS
Use Grep tool on `/claude/logs.txt` to find:
- File history: `grep "filename.py" /claude/logs.txt`
- Recent creates: `grep "\[CREATE\]" /claude/logs.txt | tail -20`
- Recent modifies: `grep "\[MODIFY\]" /claude/logs.txt | tail -20`
- Specific date: `grep "2025-06-23" /claude/logs.txt`
- Deprecated files: `grep "\[DEPRECATE\]" /claude/logs.txt`
- **View recent logs**: `tail -n 30 /claude/logs.txt`

## 📊 CONTEXT BUILDING
The logs provide:
- **Change history** for any file
- **Development timeline** 
- **Pattern detection** (repeated modifications)
- **Deprecation tracking**
- **Collaboration context** (who changed what when)

## 📝 LOG FORMAT & INSTRUCTIONS
**When logging changes:**
- Format: `[YYYY-MM-DD HH:MM] [ACTION] /path/to/file - commit message`
- Actions: CREATE, MODIFY, DELETE, RENAME (old → new), DEPRECATE
- Always APPEND to `/claude/logs.txt` (never overwrite)
- Use bash: `echo "[$(date +%Y-%m-%d' '%H:%M)] [ACTION] /path - message" >> /claude/logs.txt`

## 🎯 TASK TREE MANAGEMENT
**Project-level tasks in `/claude/todo.txt`:**
- XML format with nested task trees
- Break complex tasks into subtasks
- Status: pending | in_progress | completed | blocked
- Complete all subtasks before marking parent complete
- When solving subtasks, continue until parent task is resolved
- Update XML structure as tasks progress

### Task Tree XML Structure
Each task element can have:
- `id`: unique identifier (e.g., "1", "backend-auth", "feature.1.2")
- `status`: pending | in_progress | completed | blocked
- `priority`: high | medium | low
- `description`: what needs to be done
- `subtasks`: child tasks that must be completed first

Example format:
```xml
<task id="feature-1" status="in_progress" priority="high">
  <description>Implement user authentication system</description>
  <subtasks>
    <task id="feature-1.1" status="completed" priority="high">
      <description>Create user model with UUID</description>
    </task>
    <task id="feature-1.2" status="pending" priority="medium">
      <description>Add login/logout endpoints</description>
    </task>
  </subtasks>
</task>
```

### Task Completion Rules
1. Mark task "in_progress" when starting work
2. Complete ALL subtasks before marking parent complete
3. If blocked, add reason in description
4. When completing a task, check if parent can be completed
5. Always update todo.txt when task status changes

### Task Breakdown Guidelines
**Complex tasks should be broken down when:**
- Task requires multiple distinct operations
- Task spans multiple files or modules
- Task has dependencies on other components
- Task can be parallelized into independent subtasks

**Atomic tasks (no breakdown needed):**
- Single file edits
- Simple function implementations
- Configuration changes
- Documentation updates

---

## Project Overview

MOGRN (Microbial Opsin Generic Residue Numbering) is a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted microbial opsin structures. The project implements Generic Residue Numbering (GRN) to standardize structural comparisons across diverse opsin structures.

**Key Components:**
- Structure data processing using the Protos framework
- Multiple sequence/structure alignment with GRN assignment
- RMSD-based structural comparison and error analysis
- Transmembrane helix identification and annotation
- Visualization of structural relationships and conservation patterns

## Project Structure

```
MOGRN/
├── property/              # Property data (mo_exp.csv, helices.json)
├── structures/            # Structure files (CIF format)
│   ├── hideaki_exp/      # Experimental structures from Hideaki dataset
│   ├── hideaki_pred/     # Predicted structures from Hideaki dataset
│   └── mo_pred/          # Predicted microbial opsin structures
├── src/                  # Core modules
│   ├── data_processing.py         # Structure loading and preprocessing
│   ├── structure_comparison.py    # RMSD calculations and alignment
│   ├── helix_analysis.py          # Helix identification and annotation
│   ├── error_analysis.py          # Error analysis between structures
│   ├── msa_grn.py                # MSA with GRN assignment
│   ├── assign_grns.py            # GRN assignment workflow
│   └── visualization_functions.py # Plotting and visualization
├── protos/               # Protos framework (required dependency)
├── animation/            # Animation generation tools
└── claude/               # Claude AI session management files
```

## Key Commands

### Main Workflow
```bash
# Step 1: Initialize data infrastructure
python prepare_data_fixed.py

# Step 2: Create YAML configurations for sequences
python prepare_yaml.py

# Step 3: Run complete analysis pipeline
python opsin_analysis_workflow.py

# Step 4: Generate visualizations
python plot.py
```

### Backend Commands

```bash
# Run analysis with custom output directory
python opsin_analysis_workflow.py --output-dir custom_output

# Run without cache
python opsin_analysis_workflow.py --no-cache

# Generate plots with custom directories
python plot.py --input-dir custom_output --output-dir custom_figures

# Count chains in structures
python count_chains.py

# Debug cache functionality
python debug_cache.py
```

### Linting and Type Checking
```bash
# Python linting with flake8
flake8 src/ --max-line-length=120 --ignore=E501,W503

# Type checking with mypy
mypy src/ --ignore-missing-imports

# Run tests
pytest tests/
```

## Environment Setup

1. **Install Protos Framework** (Required):
   ```bash
   git clone https://github.com/flurinh/protos.git
   cd protos
   pip install -e .
   cd ..
   ```

2. **Install MOGRN dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Set up data directories**:
   - Create `property/` and `structures/` directories
   - Populate with required data files (see README.md)

## Important Notes

- **Path Management**: Always use absolute paths when working with file operations
- **Protos Integration**: The project heavily depends on Protos framework for structure processing
- **Memory Usage**: Large datasets may require significant memory; use caching where possible
- **File Formats**: Structures must be in CIF format; sequences in CSV/YAML
- **Chain Selection**: Most analyses focus on chain 'A' by default

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test modules
pytest tests/test_imports.py -v
pytest src/test_protos.py -v

# Test animation functionality
python animation/test_single_frame.py
python animation/test_carbon_numbering.py
```