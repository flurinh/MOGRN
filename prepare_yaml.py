import pandas as pd
import os
import re
import yaml
import matplotlib.pyplot as plt
from ruamel.yaml import YAML

# Load the dataset
df = pd.read_excel('property/mo_exp_ST1.csv', index_col=0)

# Clean sequences by removing spaces
df['seq'] = df['seq'].str.replace(' ', '')

# Reset index
df = df.reset_index(drop=True)

# Check for duplicate opsin names
duplicates = df[df.duplicated(subset=['opsin name'], keep=False)]
if not duplicates.empty:
    print("Found duplicates in 'opsin name' column:")
    print(duplicates[['opsin name']])
else:
    print("No duplicates found in 'opsin name' column")

# Check for special characters in short_name
special_char_pattern = r'[^\w]|[()]|\.'
special_chars = df[df['short_name'].notna() & df['short_name'].str.contains(special_char_pattern, regex=True)]

if not special_chars.empty:
    print("\nFound special characters in 'short_name' column:")
    for idx, row in special_chars.iterrows():
        print(f"Index: {idx}, short_name: '{row['short_name']}', opsin name: '{row['opsin name']}'")
else:
    print("\nNo special characters found in 'short_name' column")

# Fill missing short_name and display_name with opsin name
df.loc[df['short_name'].isna(), 'short_name'] = df.loc[df['short_name'].isna(), 'opsin name']
df.loc[df['display_name'].isna(), 'display_name'] = df.loc[df['display_name'].isna(), 'opsin name']

# Replace special characters with underscores in both columns
def replace_special_chars(text):
    if pd.isna(text):  # Handle NaN values
        return text
    return re.sub(special_char_pattern, '_', str(text))

df['short_name'] = df['short_name'].apply(replace_special_chars)
df['display_name'] = df['display_name'].apply(replace_special_chars)

# Print sample of updated rows
print("Sample of updated rows:")
sample_rows = df[['opsin name', 'short_name', 'display_name']].head(5)
print(sample_rows)

# Create the outputs directory if it doesn't exist
os.makedirs('yaml_configs/mo_folding6', exist_ok=True)

# Initialize YAML writer
yaml_writer = YAML()
yaml_writer.indent(mapping=2, sequence=4, offset=2)
yaml_writer.preserve_quotes = True
yaml_writer.width = 1000  # To prevent line wrapping

# Loop through each row in the dataframe to create YAML files
for index, row in df.iterrows():
    # Get the short_name for the filename
    short_name = row['short_name']
    
    # Clean the sequence by removing any whitespaces
    full_sequence = row['seq'].replace(' ', '')
    
    # Check if start and end values exist and are valid
    has_valid_range = False
    if pd.notna(row['start']) and pd.notna(row['end']):
        start = int(row['start']) if int(row['start']) >= 1 else 1
        end = int(row['end']) if int(row['end']) <= len(full_sequence) else len(full_sequence)
        if start < end:
            sequence = full_sequence[start-1:end]  # Adjust for 0-based indexing
            has_valid_range = True
    
    # If no valid range, use full sequence
    if not has_valid_range:
        sequence = full_sequence
        if len(sequence) > 400:
            print(f"WARNING: Long sequence without range - {short_name}: length {len(sequence)}")
    
    # Create the YAML structure
    yaml_data = {
        'version': 1,
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': sequence
                }
            },
            {
                'ligand': {
                    'id': 'B',
                    'smiles': 'CC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C'
                }
            }
        ]
    }
    
    # Create filename based on short_name
    filename = f'yaml_configs/mo_folding6/{short_name}.yaml'
    
    # Write the file with proper formatting
    with open(filename, 'w') as file:
        yaml_writer.dump(yaml_data, file)
    
    # Add the comment after "version: 1"
    with open(filename, 'r') as file:
        content = file.read()
    
    content = content.replace('version: 1', 'version: 1  # Optional, defaults to 1')
    
    with open(filename, 'w') as file:
        file.write(content)
    
    print(f"Created YAML file: {filename}")

print(f"Finished creating {len(df)} YAML files in yaml_configs/mo_folding6/")

# Create special YAML file with two ligands for Tara_RRB
tara_row = df[df['short_name'] == 'Tara_RRB']

if not tara_row.empty:
    # Get the sequence and clean it
    full_sequence = tara_row['seq'].iloc[0].replace(' ', '')
    
    # Check if start and end values exist and are valid
    has_valid_range = False
    if pd.notna(tara_row['start'].iloc[0]) and pd.notna(tara_row['end'].iloc[0]):
        start = int(tara_row['start'].iloc[0]) if int(tara_row['start'].iloc[0]) >= 1 else 1
        end = int(tara_row['end'].iloc[0]) if int(tara_row['end'].iloc[0]) <= len(full_sequence) else len(full_sequence)
        if start < end:
            sequence = full_sequence[start-1:end]  # Adjust for 0-based indexing
            has_valid_range = True
    
    # If no valid range, use full sequence
    if not has_valid_range:
        sequence = full_sequence
        if len(sequence) > 400:
            print(f"WARNING: Long sequence without range in special file - Tara_RRB: length {len(sequence)}")
    
    # Create the special YAML structure with two ligands
    yaml_data = {
        'version': 1,
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': sequence
                }
            },
            {
                'ligand': {
                    'id': 'B',
                    'smiles': 'CC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C'
                }
            },
            {
                'ligand': {
                    'id': 'C',
                    'smiles': 'CC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C'
                }
            }
        ]
    }
    
    # Create filename
    filename = 'yaml_configs/mo_folding6/Tara_RRB_2_retinals.yaml'
    
    # Write the file with proper formatting
    with open(filename, 'w') as file:
        yaml_writer.dump(yaml_data, file)
    
    # Add the comment after "version: 1"
    with open(filename, 'r') as file:
        content = file.read()
    
    content = content.replace('version: 1', 'version: 1  # Optional, defaults to 1')
    
    with open(filename, 'w') as file:
        file.write(content)
    
    print(f"Created special YAML file with two ligands: {filename}")
else:
    print("Error: Could not find 'Tara_RRB' in the dataframe")

# Validate amino acid sequences
def validate_amino_acid_sequence(sequence):
    # Define valid amino acid characters (standard 20 amino acids)
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    # Remove any whitespaces first
    sequence = sequence.replace(' ', '')
    
    # Check if any invalid characters exist
    invalid_chars = set(sequence) - valid_amino_acids
    
    if invalid_chars:
        return False, invalid_chars
    else:
        return True, set()

# Check all sequences in the dataframe
invalid_sequences = []

for index, row in df.iterrows():
    sequence = row['seq']
    short_name = row['short_name']
    
    is_valid, invalid_chars = validate_amino_acid_sequence(sequence)
    
    if not is_valid:
        invalid_sequences.append({
            'index': index,
            'short_name': short_name,
            'invalid_chars': ''.join(invalid_chars)
        })

# Display results of validation
if invalid_sequences:
    print(f"Found {len(invalid_sequences)} sequences with invalid amino acid characters:")
    for seq in invalid_sequences:
        print(f"Row {seq['index']}, short_name: '{seq['short_name']}', Invalid characters: '{seq['invalid_chars']}'")
else:
    print("All sequences contain valid amino acid characters only.")

# Calculate and analyze sequence lengths
df['len'] = df['end'] - df['start']
df['len'] = df['len'].fillna(0)
df['seq_len'] = df['seq'].str.replace(' ', '').str.len()

# Add a new column for the length of the extracted sequence segment
df['extracted_seq_len'] = None
for idx, row in df.iterrows():
    full_sequence = row['seq'].replace(' ', '')
    if pd.notna(row['start']) and pd.notna(row['end']):
        start = int(row['start']) if int(row['start']) >= 1 else 1
        end = int(row['end']) if int(row['end']) <= len(full_sequence) else len(full_sequence)
        if start < end:
            df.at[idx, 'extracted_seq_len'] = end - start + 1
    else:
        df.at[idx, 'extracted_seq_len'] = len(full_sequence)

# Find entries with length > 400
long_entries = df[df['seq_len'] > 400]

# Print short names of these entries along with both lengths
if not long_entries.empty:
    print(f"\nFound {len(long_entries)} entries with length > 400:")
    print("FORMAT: short_name: calculated_length (end-start), actual_sequence_length, extracted_sequence_length")
    for _, row in long_entries.iterrows():
        print(f"{row['short_name']}: {row['len']}, {row['seq_len']}, {row['extracted_seq_len']}")
else:
    print("No entries found with length > 400")

# Summarize which entries are using full sequence vs extracted segment
with_range = df[pd.notna(df['start']) & pd.notna(df['end'])]
print(f"\nEntries with defined range (start:end): {len(with_range)} of {len(df)}")
print(f"Entries using full sequence: {len(df) - len(with_range)} of {len(df)}")

# Create histograms of sequence lengths
plt.figure(figsize=(12, 8))

# First subplot for full sequences
plt.subplot(2, 1, 1)
plt.hist(df['seq_len'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
plt.title('Distribution of Full Sequence Lengths')
plt.xlabel('Full Sequence Length')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# Add stats for full sequences
mean_len = df['seq_len'].mean()
median_len = df['seq_len'].median()
min_len = df['seq_len'].min()
max_len = df['seq_len'].max()

stats_text = f"Mean: {mean_len:.1f}\nMedian: {median_len:.1f}\nMin: {min_len}\nMax: {max_len}"
plt.annotate(stats_text, xy=(0.75, 0.75), xycoords='axes fraction', 
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

# Second subplot for extracted sequences
plt.subplot(2, 1, 2)
plt.hist(df['extracted_seq_len'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
plt.title('Distribution of Extracted Sequence Lengths (Using start:end when available)')
plt.xlabel('Extracted Sequence Length')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# Add stats for extracted sequences
mean_extracted = df['extracted_seq_len'].mean()
median_extracted = df['extracted_seq_len'].median()
min_extracted = df['extracted_seq_len'].min()
max_extracted = df['extracted_seq_len'].max()

stats_text = f"Mean: {mean_extracted:.1f}\nMedian: {median_extracted:.1f}\nMin: {min_extracted}\nMax: {max_extracted}"
plt.annotate(stats_text, xy=(0.75, 0.75), xycoords='axes fraction', 
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

plt.tight_layout()

plt.savefig('sequence_length_histogram.png')
print("Histogram of sequence lengths has been created and saved as 'sequence_length_histogram.png'")

# Save the processed dataframe
df.to_csv('property/mo_exp.csv')
print("Saved processed dataframe to 'property/mo_exp.csv'")