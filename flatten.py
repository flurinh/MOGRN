import os
import shutil
import re

OUTPUT_ROOT = "output2"
DEST_DIR = "flat_output2"


def to_long_path(path: str) -> str:
    """Convert a path to an extended-length path on Windows."""
    if os.name == "nt":
        abs_path = os.path.abspath(path)
        if not abs_path.startswith("\\\\?\\"):
            return "\\\\?\\" + abs_path
        return abs_path
    return path


def clean_sample_name(folder_name):
    """
    Remove the trailing '.yaml' extension (if present) and any unwanted substrings like '_s199f'
    from the folder name.
    """
    base, ext = os.path.splitext(folder_name)
    sample = base if ext.lower() == ".yaml" else folder_name
    return sample


def construct_cif_path(clean_folder_path, sample_name):
    """
    Construct the expected CIF file path using the cleaned folder path and sample name.

    Expected structure:
      <clean_folder_path>/boltz_results_<sample_name>/predictions/<sample_name>/<sample_name>_model_0.cif
    """
    cif_rel_path = os.path.join(
        f"boltz_results_{sample_name}",
        "predictions",
        sample_name,
        f"{sample_name}_model_0.cif"
    )
    return os.path.join(clean_folder_path, cif_rel_path)


def main():
    os.makedirs(DEST_DIR, exist_ok=True)

    for entry in os.listdir(OUTPUT_ROOT):
        subfolder_path = os.path.join(OUTPUT_ROOT, entry)
        if not os.path.isdir(subfolder_path):
            continue

        # Clean the folder name and folder path.
        sample_name = clean_sample_name(entry)
        # Also clean the actual folder path by removing '_s199f'
        clean_folder_path = subfolder_path

        # Construct the expected CIF file path.
        expected_cif = construct_cif_path(clean_folder_path, sample_name)

        print(f"Checking sample '{sample_name}' in folder:\n  {subfolder_path}")
        print(f"Using cleaned folder path:\n  {clean_folder_path}")
        print(f"Constructed CIF path:\n  {expected_cif}\n")

        # Check existence using the long path helper
        if os.path.exists(to_long_path(expected_cif)):
            dest_filename = f"{sample_name}_model_0.cif"
            dest_path = os.path.join(DEST_DIR, dest_filename)
            try:
                shutil.copyfile(to_long_path(expected_cif), to_long_path(dest_path))
                print(f"Copied: {expected_cif} -> {dest_path}\n")
            except Exception as e:
                print(f"Error copying {expected_cif}: {e}\n")
        else:
            print(f"CIF file NOT found for sample '{sample_name}' in folder:\n  {subfolder_path}\n")


if __name__ == "__main__":
    main()
