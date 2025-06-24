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


print("Total number of structures:", len(processed_structures))
print(processed_structures.keys())
print(list(processed_structures['4PXK'].keys()))
print(processed_structures['4PXK']['structure_type'])
print(processed_structures['4PXK']['base_name'])
print(processed_structures['HsHR_model_0']['base_name'])
print(processed_structures['R2ACR_J315_refine8']['base_name'])

property_csv_path = project_dir / 'property' / 'mo_exp.csv'
if property_csv_path.exists():
    property_data = load_opsin_property_data(property_csv_path, processed_structures)

print("number of property entries:", len(property_data['properties']))
original = property_data["structure_mapping"]
inverted = {v: k for k, v in original.items()} # inverted mapping (keys are predictions)



print(property_data['properties']['3UG9'])


def load_grn_tables_data():  # Matching notebook name
    grn_tables_pkl = GRN_TABLES_DIR / 'grn_tables.pkl'  # Corrected path
    return load_cached_data(grn_tables_pkl, "GRN tables data")


grn_tables = load_grn_tables_data()
print(grn_tables.keys())
print(grn_tables['residue_table'].head(2))
print(grn_tables['helix_pivot_columns'])

# Count properties for experimental and predicted structures
import matplotlib.pyplot as plt
from collections import Counter
from src.opsin_color_scheme import get_categorical_colors, OPSIN_COLORS

experimental_properties = {}
predicted_properties = {}

for struct_id, properties in property_data['properties'].items():
    # Check if structure exists in processed_structures
    struct_id = struct_id.replace("_smile", "")
    if (struct_id in processed_structures) & (struct_id not in inverted.keys()):
        struct_type = processed_structures[struct_id]['structure_type']
        if 'model' in processed_structures[struct_id]['base_name']:
            struct_type = 'predicted'
        if struct_type == 'experimental':
            experimental_properties[struct_id] = properties
        elif struct_type == 'predicted':
            predicted_properties[struct_id] = properties
    else:
        print("not in structures", struct_id)

print(f"\nTotal experimental structures with properties: {len(experimental_properties)}")
print(f"Total predicted structures with properties: {len(predicted_properties)}")

# Get all unique property keys
all_property_keys = set()
for props in experimental_properties.values():
    all_property_keys.update(props.keys())
for props in predicted_properties.values():
    all_property_keys.update(props.keys())

print(f"\nUnique property types found: {sorted(all_property_keys)}")

# Count occurrences of each property value for both experimental and predicted
property_counts = {}

for prop_key in sorted(all_property_keys):
    # Count for experimental
    exp_values = []
    for props in experimental_properties.values():
        if prop_key in props:
            exp_values.append(props[prop_key])
    
    # Count for predicted
    pred_values = []
    for props in predicted_properties.values():
        if prop_key in props:
            pred_values.append(props[prop_key])
    
    property_counts[prop_key] = {
        'experimental': Counter(exp_values),
        'predicted': Counter(pred_values)
    }

# Create bar plots for each property
fig_dir = PROJECT_DIR / 'opsin_output' / 'figures' / 'property_analysis'
fig_dir.mkdir(parents=True, exist_ok=True)

for prop_key, counts in property_counts.items():
    # Get all unique values across both experimental and predicted
    all_values = set()
    exp_counts = counts['experimental']
    pred_counts = counts['predicted']
    all_values.update(exp_counts.keys())
    all_values.update(pred_counts.keys())
    all_values = sorted(list(all_values))  # Sort for consistent ordering
    
    # Determine property type and get colors
    if prop_key.lower() == 'molecular_function' or 'function' in prop_key.lower():
        property_type = 'property1'
    elif prop_key.lower() == 'domain' or 'kingdom' in prop_key.lower():
        property_type = 'property2'
    else:
        property_type = 'property1'  # Default to warm colors
    
    # Get colors for all values
    value_colors = get_categorical_colors(all_values, property_type=property_type)
    bar_colors = [value_colors[val] for val in all_values]
    
    # Create figure with subplots for experimental and predicted
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot experimental data with consistent ordering
    if exp_counts:
        counts_list = [exp_counts.get(val, 0) for val in all_values]
        
        bars1 = ax1.bar(range(len(all_values)), counts_list, color=bar_colors)
        ax1.set_xticks(range(len(all_values)))
        ax1.set_xticklabels(all_values, rotation=45, ha='right')
        ax1.set_xlabel(prop_key)
        ax1.set_ylabel('Count')
        ax1.set_title(f'Experimental Structures - {prop_key}')
        ax1.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar, count in zip(bars1, counts_list):
            if count > 0:  # Only show label if count > 0
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{count}', ha='center', va='bottom', fontsize=9)
    
    # Plot predicted data with consistent ordering
    if pred_counts:
        counts_list = [pred_counts.get(val, 0) for val in all_values]
        
        bars2 = ax2.bar(range(len(all_values)), counts_list, color=bar_colors)
        ax2.set_xticks(range(len(all_values)))
        ax2.set_xticklabels(all_values, rotation=45, ha='right')
        ax2.set_xlabel(prop_key)
        ax2.set_ylabel('Count')
        ax2.set_title(f'Predicted Structures - {prop_key}')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar, count in zip(bars2, counts_list):
            if count > 0:  # Only show label if count > 0
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{count}', ha='center', va='bottom', fontsize=9)
    
    plt.suptitle(f'Property Distribution: {prop_key}', fontsize=16)
    plt.tight_layout()
    
    # Save figure
    filename = f'property_{prop_key.replace(" ", "_").replace("/", "_")}_distribution.png'
    filepath = fig_dir / filename
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {filepath}")
    plt.close()

# Create a combined summary plot showing counts of structures with each property
fig, ax = plt.subplots(figsize=(10, 6))

prop_names = sorted(all_property_keys)
exp_has_prop = [sum(1 for props in experimental_properties.values() if prop in props) for prop in prop_names]
pred_has_prop = [sum(1 for props in predicted_properties.values() if prop in props) for prop in prop_names]

x = range(len(prop_names))
width = 0.35

bars_exp = ax.bar([i - width/2 for i in x], exp_has_prop, width, label='Experimental', color='steelblue')
bars_pred = ax.bar([i + width/2 for i in x], pred_has_prop, width, label='Predicted', color='coral')

# Add value labels on bars
for bar, count in zip(bars_exp, exp_has_prop):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{count}', ha='center', va='bottom', fontsize=9)

for bar, count in zip(bars_pred, pred_has_prop):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{count}', ha='center', va='bottom', fontsize=9)

ax.set_xlabel('Property Type')
ax.set_ylabel('Number of Structures')
ax.set_title('Number of Structures with Each Property Type')
ax.set_xticks(x)
ax.set_xticklabels(prop_names, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
summary_filepath = fig_dir / 'property_coverage_summary.png'
plt.savefig(summary_filepath, dpi=300, bbox_inches='tight')
print(f"Saved summary plot: {summary_filepath}")
plt.close()

print("\nProperty analysis complete!")

# Create combined total counts plot for each property
print("\nCreating combined total counts plots...")

for prop_key, counts in property_counts.items():
    # Combine experimental and predicted counts
    exp_counts = counts['experimental']
    pred_counts = counts['predicted']
    
    # Get all unique values and their total counts
    total_counts = Counter()
    for val, count in exp_counts.items():
        total_counts[val] += count
    for val, count in pred_counts.items():
        total_counts[val] += count
    
    if total_counts:
        # Sort by value for consistent ordering
        sorted_items = sorted(total_counts.items())
        values = [item[0] for item in sorted_items]
        counts_list = [item[1] for item in sorted_items]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Use same colors as in individual plots
        if prop_key.lower() == 'molecular_function' or 'function' in prop_key.lower():
            property_type = 'property1'
        elif prop_key.lower() == 'domain' or 'kingdom' in prop_key.lower():
            property_type = 'property2'
        else:
            property_type = 'property1'
        
        value_colors = get_categorical_colors(values, property_type=property_type)
        bar_colors = [value_colors[val] for val in values]
        
        bars = ax.bar(range(len(values)), counts_list, color=bar_colors)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(values, rotation=45, ha='right')
        ax.set_xlabel(prop_key)
        ax.set_ylabel('Total Count')
        ax.set_title(f'Total Structures (Experimental + Predicted) - {prop_key}')
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar, count in zip(bars, counts_list):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        
        # Save figure
        filename = f'property_{prop_key.replace(" ", "_").replace("/", "_")}_total.png'
        filepath = fig_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved total plot: {filepath}")
        plt.close()

print("\nAll plots generated successfully!")

# Create combination heatmap for molecular function vs domain
print("\nCreating property combination heatmap...")

# Combine all properties from experimental and predicted
all_properties = {}
for struct_id, props in experimental_properties.items():
    all_properties[struct_id] = props
for struct_id, props in predicted_properties.items():
    if struct_id not in all_properties:  # Avoid duplicates
        all_properties[struct_id] = props

# Find molecular function and domain properties
function_key = None
domain_key = None

for key in all_property_keys:
    if 'molecular_function' in key.lower() or 'function' in key.lower():
        function_key = key
    elif 'domain' in key.lower() or 'kingdom' in key.lower():
        domain_key = key

if function_key and domain_key:
    print(f"Found function key: {function_key}")
    print(f"Found domain key: {domain_key}")
    
    # Create combination counts
    import pandas as pd
    import seaborn as sns
    
    combinations = []
    # Also create filtered combinations for natural domains and known functions
    natural_combinations = []
    
    for struct_id, props in all_properties.items():
        if function_key in props and domain_key in props:
            combinations.append({
                'function': props[function_key],
                'domain': props[domain_key]
            })
            # Filter for natural domains (exclude Synthetic) and known functions (exclude Unknown)
            if props[domain_key].lower() != 'synthetic' and props[function_key].lower() != 'unknown':
                natural_combinations.append({
                    'function': props[function_key],
                    'domain': props[domain_key]
                })
    
    # Create DataFrame and pivot table for all combinations
    df_combinations = pd.DataFrame(combinations)
    pivot_table = pd.crosstab(df_combinations['function'], df_combinations['domain'])
    
    # Create DataFrame and pivot table for natural combinations only
    df_natural_combinations = pd.DataFrame(natural_combinations)
    if len(df_natural_combinations) > 0:
        natural_pivot_table = pd.crosstab(df_natural_combinations['function'], df_natural_combinations['domain'])
    else:
        natural_pivot_table = pd.DataFrame()
    
    # Sort rows and columns
    pivot_table = pivot_table.sort_index()
    pivot_table = pivot_table[sorted(pivot_table.columns)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create custom annotations - only show non-zero values
    annot_data = pivot_table.astype(str)
    annot_data[pivot_table == 0] = ''  # Empty string for zero values
    
    # Use grayscale colormap from opsin_color_scheme
    from src.opsin_color_scheme import RMSD_WHITE_TO_DARKGRAY_CMAP
    
    # Create heatmap
    sns.heatmap(pivot_table, 
                annot=annot_data,  # Use custom annotations
                fmt='',     # Use string format (no additional formatting)
                cmap=RMSD_WHITE_TO_DARKGRAY_CMAP,  # White to dark gray colormap
                linewidths=0.5,
                linecolor='gray',
                cbar_kws={'label': 'Count'},
                vmin=0,  # Ensure 0 is white
                ax=ax)
    
    ax.set_xlabel(domain_key, fontsize=12)
    ax.set_ylabel(function_key, fontsize=12)
    ax.set_title('Combination of Molecular Function and Domain Properties', fontsize=14)
    
    # Rotate labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    # Save figure
    heatmap_filepath = fig_dir / 'property_combination_heatmap.png'
    plt.savefig(heatmap_filepath, dpi=300, bbox_inches='tight')
    print(f"Saved combination heatmap: {heatmap_filepath}")
    plt.close()
    
    # Also create a summary statistics
    print("\nCombination statistics:")
    print(f"Total structures with both properties: {len(combinations)}")
    print(f"Unique functions: {pivot_table.shape[0]}")
    print(f"Unique domains: {pivot_table.shape[1]}")
    print(f"\nTop 5 most common combinations:")
    flat_combinations = [(f, d, pivot_table.loc[f, d]) 
                         for f in pivot_table.index 
                         for d in pivot_table.columns 
                         if pivot_table.loc[f, d] > 0]
    flat_combinations.sort(key=lambda x: x[2], reverse=True)
    for i, (func, dom, count) in enumerate(flat_combinations[:5]):
        print(f"  {i+1}. {func} + {dom}: {count} structures")
    
    # Read existing diversity.md to preserve literature review section
    diversity_file = PROJECT_DIR / 'diversity.md'
    existing_content = ""
    preserve_from_marker = "## Comparison with Published Literature"
    
    if diversity_file.exists():
        with open(diversity_file, 'r', encoding='utf-8') as f:
            content = f.read()
            if preserve_from_marker in content:
                # Find the start of the section to preserve
                preserve_index = content.find(preserve_from_marker)
                existing_content = content[preserve_index:]
    
    # Write diversity analysis to file
    with open(diversity_file, 'w', encoding='utf-8') as f:
        f.write("# Microbial Opsin Diversity Analysis\n\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Dataset Overview\n\n")
        f.write(f"- Total structures analyzed: {len(all_properties)}\n")
        f.write(f"- Experimental structures: {len(experimental_properties)}\n")
        f.write(f"- Predicted structures: {len(predicted_properties)}\n")
        f.write(f"- Structures with both function and domain annotations: {len(combinations)}\n\n")
        
        f.write("## Molecular Function Distribution\n\n")
        function_counts = df_combinations['function'].value_counts()
        for func, count in function_counts.items():
            percentage = (count / len(combinations)) * 100
            f.write(f"- **{func}**: {count} structures ({percentage:.1f}%)\n")
        
        f.write("\n## Domain Distribution\n\n")
        domain_counts = df_combinations['domain'].value_counts()
        for dom, count in domain_counts.items():
            percentage = (count / len(combinations)) * 100
            f.write(f"- **{dom}**: {count} structures ({percentage:.1f}%)\n")
        
        f.write("\n## Function-Domain Combinations\n\n")
        f.write("### Top 20 Most Common Combinations\n\n")
        f.write("| Rank | Molecular Function | Domain | Count | Percentage |\n")
        f.write("|------|-------------------|---------|--------|------------|\n")
        for i, (func, dom, count) in enumerate(flat_combinations[:20]):
            percentage = (count / len(combinations)) * 100
            f.write(f"| {i+1} | {func} | {dom} | {count} | {percentage:.1f}% |\n")
        
        f.write("\n### Complete Combination Matrix\n\n")
        f.write("Rows: Molecular Function, Columns: Domain\n\n")
        
        # Write the pivot table as markdown
        f.write("| Function / Domain |")
        for col in pivot_table.columns:
            f.write(f" {col} |")
        f.write("\n")
        
        f.write("|-------------------|")
        for _ in pivot_table.columns:
            f.write("--------|")
        f.write("\n")
        
        for idx in pivot_table.index:
            f.write(f"| {idx} |")
            for col in pivot_table.columns:
                val = pivot_table.loc[idx, col]
                if val > 0:
                    f.write(f" {val} |")
                else:
                    f.write(" - |")
            f.write("\n")
        
        f.write("\n## Diversity Metrics\n\n")
        # Calculate Shannon diversity index for functions
        from scipy.stats import entropy
        function_probs = function_counts / function_counts.sum()
        function_shannon = entropy(function_probs)
        
        # Calculate Shannon diversity index for domains
        domain_probs = domain_counts / domain_counts.sum()
        domain_shannon = entropy(domain_probs)
        
        f.write(f"- **Shannon Diversity Index (Molecular Function)**: {function_shannon:.3f}\n")
        f.write(f"- **Shannon Diversity Index (Domain)**: {domain_shannon:.3f}\n")
        f.write(f"- **Number of unique function-domain combinations**: {(pivot_table > 0).sum().sum()}\n")
        f.write(f"- **Theoretical maximum combinations**: {len(function_counts) * len(domain_counts)}\n")
        f.write(f"- **Combination coverage**: {((pivot_table > 0).sum().sum() / (len(function_counts) * len(domain_counts)) * 100):.1f}%\n")
        
        # Add natural combinations analysis
        if len(natural_pivot_table) > 0:
            # Count natural domains and functions
            natural_function_counts = df_natural_combinations['function'].value_counts()
            natural_domain_counts = df_natural_combinations['domain'].value_counts()
            
            f.write("\n### Natural Combinations Only (Excluding Synthetic Domain and Unknown Function)\n\n")
            f.write(f"- **Natural functions**: {len(natural_function_counts)} (excluding Unknown)\n")
            f.write(f"- **Natural domains**: {len(natural_domain_counts)} (excluding Synthetic)\n")
            f.write(f"- **Natural function-domain combinations observed**: {(natural_pivot_table > 0).sum().sum()}\n")
            f.write(f"- **Theoretical natural combinations**: {len(natural_function_counts) * len(natural_domain_counts)} (6 functions × 4 domains)\n")
            f.write(f"- **Natural combination coverage**: {((natural_pivot_table > 0).sum().sum() / (len(natural_function_counts) * len(natural_domain_counts)) * 100):.1f}%\n")
        
        # Append the preserved content (literature review)
        if existing_content:
            f.write("\n")
            f.write(existing_content)
        
    print(f"\nDiversity analysis written to: {diversity_file}")
    print(f"Natural combinations: {len(natural_combinations)} structures (excluding Synthetic domain and Unknown function)")
    if len(natural_pivot_table) > 0:
        print(f"Natural coverage: {(natural_pivot_table > 0).sum().sum()} / {6 * 4} = {((natural_pivot_table > 0).sum().sum() / 24 * 100):.1f}%")
else:
    print("Could not find both molecular function and domain properties in the data.")

print("\nProperty combination analysis complete!")

