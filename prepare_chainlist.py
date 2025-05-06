import os

struct_dir = "all_structs/"
output_file = "chainlist.txt"
suffix = ".cif"

base_filenames = []
try:
    all_entries = os.listdir(struct_dir)
except FileNotFoundError:
    print(f"Error: Directory not found: {struct_dir}")
    exit() # Or raise an exception

for entry_name in all_entries:
    full_path = os.path.join(struct_dir, entry_name) # Needed for isfile check
    if os.path.isfile(full_path) and entry_name.endswith(suffix):
        base_name, ext = os.path.splitext(entry_name)
        if ext == suffix: # Double check extension
            base_filenames.append(base_name)

# Write to file
try:
    with open(output_file, 'w') as f:
        for name in base_filenames:
            f.write(name + "\n")
    print(f"Successfully created {output_file} with {len(base_filenames)} entries.")
except IOError:
    print(f"Error: Could not write to file {output_file}")

# then we run in Ubuntu: USalign -dir all_structs/ chainlist.txt -mol prot -outfmt 2