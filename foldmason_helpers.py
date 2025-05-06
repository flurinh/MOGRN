import subprocess
import shlex
import os
from tqdm import tqdm

def run_command(command, use_wsl=True):
    """
    Execute a shell command and return its stdout output.
    If use_wsl is True, the command is run via WSL in an interactive login shell,
    ensuring that your PATH (and other settings) are properly loaded.

    Raises subprocess.CalledProcessError if the command fails.
    """
    if use_wsl:
        # Use bash -lic to force an interactive login shell in WSL.
        command = f"wsl bash -lic {shlex.quote(command)}"
    try:
        result = subprocess.run(shlex.split(command),
                                capture_output=True,
                                text=True,
                                check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {command}")
        print(f"Error: {e.stderr.strip()}")
        raise


def run_easy_msa(input_files, output_prefix, tmp_folder, report_mode=1, precluster=False):
    """
    Run the 'easy-msa' command.

    Usage:
      foldmason easy-msa <PDB/mmCIF files> <output_prefix> <tmp_folder> [--report-mode <mode>] [--precluster]

    Parameters:
      input_files (list of str): List of input PDB/mmCIF file paths.
      output_prefix (str): Output file prefix.
      tmp_folder (str): Path to a temporary folder.
      report_mode (int, optional): Report mode (default is 1).
      precluster (bool, optional): If True, add the '--precluster' flag.

    Returns:
      str: The stdout output of the command.
    """
    files_str = " ".join(input_files)
    cmd = f"foldmason easy-msa {files_str} {output_prefix} {tmp_folder} --report-mode {report_mode}"
    cmd = f"foldmason easy-msa {files_str} {output_prefix} {tmp_folder} --report-mode {report_mode}"

    if precluster:
        cmd += " --precluster"
    return run_command(cmd)


def run_convertalis(input_db, output_format, output_file, extra_args=""):
    """
    Run the 'convertalis' command.

    Usage:
      foldmason convertalis <input_db> <output_format> <output_file> [<extra_args>]

    Parameters:
      input_db (str): Path to the alignment DB.
      output_format (str): Desired output format (e.g., 'BLAST-tab', 'SAM').
      output_file (str): Path for the converted output.
      extra_args (str, optional): Any extra command-line arguments.

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason convertalis {input_db} {output_format} {output_file} {extra_args}".strip()
    return run_command(cmd)


def run_structuremsa(structure_db, output_prefix):
    """
    Run the 'structuremsa' command.

    Usage:
      foldmason structuremsa <structure_db> <output_prefix>

    Parameters:
      structure_db (str): Path to the structure database.
      output_prefix (str): Output file prefix.

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason structuremsa {structure_db} {output_prefix}"
    return run_command(cmd)


def run_structuremsacluster(structure_db, output_prefix, cluster_threshold=0.5):
    """
    Run the 'structuremsacluster' command.

    Usage:
      foldmason structuremsacluster <structure_db> <output_prefix> [--cluster-threshold <value>]

    Parameters:
      structure_db (str): Path to the structure database.
      output_prefix (str): Output file prefix.
      cluster_threshold (float, optional): Clustering threshold (default 0.5).

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason structuremsacluster {structure_db} {output_prefix} --cluster-threshold {cluster_threshold}"
    return run_command(cmd)


def run_msa2lddt(structure_db, input_fasta):
    """
    Run the 'msa2lddt' command to calculate the LDDT score of an MSA.

    Usage:
      foldmason msa2lddt <structure_db> <input_fasta>

    Parameters:
      structure_db (str): Path to the structure database.
      input_fasta (str): Path to the input FASTA alignment file.

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason msa2lddt {structure_db} {input_fasta}"
    return run_command(cmd)


def run_msa2lddtreport(structure_db, input_fasta, output_html, guide_tree=None):
    """
    Run the 'msa2lddtreport' command to calculate the LDDT and generate an HTML report.

    Usage:
      foldmason msa2lddtreport <structure_db> <input_fasta> <output_html> [--guide-tree <guide_tree>]

    Parameters:
      structure_db (str): Path to the structure database.
      input_fasta (str): Path to the input FASTA alignment file.
      output_html (str): Path for the HTML report.
      guide_tree (str, optional): Path to the guide tree file.

    Returns:
      str: The stdout output of the command.
    """
    if guide_tree:
        cmd = f"foldmason msa2lddtreport {structure_db} {input_fasta} {output_html} --guide-tree {guide_tree}"
    else:
        cmd = f"foldmason msa2lddtreport {structure_db} {input_fasta} {output_html}"
    return run_command(cmd)


def run_msa2lddtjson(structure_db, input_fasta, output_json):
    """
    Run the 'msa2lddtjson' command to calculate LDDT and generate a JSON report.

    Usage:
      foldmason msa2lddtjson <structure_db> <input_fasta> <output_json>

    Parameters:
      structure_db (str): Path to the structure database.
      input_fasta (str): Path to the input FASTA alignment file.
      output_json (str): Path for the JSON report.

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason msa2lddtjson {structure_db} {input_fasta} {output_json}"
    return run_command(cmd)


def run_refinemsa(structure_db, input_fasta, output_fasta, refine_iters=1000):
    """
    Run the 'refinemsa' command to iteratively refine an MSA.

    Usage:
      foldmason refinemsa <structure_db> <input_fasta> <output_fasta> --refine-iters <iterations>

    Parameters:
      structure_db (str): Path to the structure database.
      input_fasta (str): Path to the input FASTA alignment file.
      output_fasta (str): Path where the refined alignment will be saved.
      refine_iters (int, optional): Number of refinement iterations (default: 1000).

    Returns:
      str: The stdout output of the command.
    """
    cmd = f"foldmason refinemsa {structure_db} {input_fasta} {output_fasta} --refine-iters {refine_iters}"
    return run_command(cmd)


if __name__ == "__main__":
    # Example usage: print the help output for FoldMason
    try:
        output = run_command("foldmason --help")
        print("FoldMason help:")
        print(output)
    except Exception as err:
        print("Error running FoldMason:", err)


def prepare_structures_for_foldmason(processed_structures_complete, tmp_dir='tmp_foldmason'):
    """
    Export processed structures to temporary CIF files for FoldMason using direct file operations
    """
    print("[INFO] Preparing structures for FoldMason alignment...")
    os.makedirs(tmp_dir, exist_ok=True)
    structure_files = []

    # Create a temporary CifProcessor for writing files
    from protos.processing.structure.struct_base_processor import CifBaseProcessor
    data_dir = os.environ.get("PROTOS_DATA_ROOT", os.path.abspath("opsin_output"))
    cp_temp = CifBaseProcessor(name="tmp_foldmason", data_root=data_dir)

    for pdb_id, data in tqdm(processed_structures_complete.items(), desc="Creating temporary CIF files"):
        if 'df_norm' not in data or data['df_norm'].empty:
            print(f"[WARNING] {pdb_id}: No df_norm available, skipping.")
            continue

        # Create a temporary dataframe with the oriented, processed structure
        temp_df = data['df_norm'].copy()

        # Ensure all necessary columns are present for CIF output
        if 'pdb_id' not in temp_df.columns:
            temp_df['pdb_id'] = pdb_id

        # Ensure atom names are properly set if missing
        required_columns = ['group', 'type_symbol', 'auth_asym_id', 'auth_atom_id']
        for col in required_columns:
            if col not in temp_df.columns:
                if col == 'group':
                    temp_df[col] = 'ATOM'
                elif col == 'type_symbol':
                    # Derive from res_atom_name if available
                    if 'res_atom_name' in temp_df.columns:
                        temp_df[col] = temp_df['res_atom_name'].str[0]
                    else:
                        temp_df[col] = 'C'  # Default to carbon
                elif col == 'auth_asym_id':
                    if 'auth_chain_id' in temp_df.columns:
                        temp_df[col] = temp_df['auth_chain_id']
                    else:
                        temp_df[col] = 'A'
                elif col == 'auth_atom_id':
                    if 'res_atom_name' in temp_df.columns:
                        temp_df[col] = temp_df['res_atom_name']
                    else:
                        temp_df[col] = 'CA'  # Default to CA atoms

        # Write directly to a temporary CIF file
        temp_file = os.path.join(tmp_dir, f"{pdb_id}.cif")

        try:
            # Use cif_handler directly for more reliable output
            from protos.io.cif_handler import write_cif_file
            write_cif_file(temp_df, temp_file)
            structure_files.append(temp_file)
        except Exception as e:
            print(f"[ERROR] Failed to create temporary CIF for {pdb_id}: {str(e)}")

    print(f"[INFO] Created {len(structure_files)} temporary CIF files for FoldMason.")
    return structure_files


def align_with_foldmason(data_dict, output_dir='output', tmp_dir='tmp_foldmason', visualize=True):
    """
    Use FoldMason for multiple structure alignment instead of pairwise alignment
    """
    print("\n=== Step 6: FoldMason Multiple Structure Alignment & GRN Assignment ===")

    processed_structures_complete = data_dict['processed_structures']

    # Make sure output directories exist
    os.makedirs(output_dir, exist_ok=True)
    tmp_path = os.path.join(output_dir, tmp_dir)
    os.makedirs(tmp_path, exist_ok=True)

    # 1. Export structures to temporary files
    structure_files = prepare_structures_for_foldmason(processed_structures_complete, tmp_path)

    if not structure_files:
        print("[ERROR] No structure files were created. Cannot proceed with FoldMason alignment.")
        return data_dict

    try:
        # 2. Run FoldMason's easy-msa for initial alignment
        output_prefix = os.path.join(output_dir, "opsin_msa")
        print(f"[INFO] Running FoldMason easy-msa on {len(structure_files)} structures...")
        run_easy_msa(structure_files, output_prefix, tmp_path, report_mode=2)

        # 3. Generate refined structural MSA
        structure_db = f"{output_prefix}_strucDB"
        msa_prefix = os.path.join(output_dir, "opsin_structural_msa")

        print("[INFO] Running structural MSA with FoldMason...")
        run_structuremsa(structure_db, msa_prefix)

        # 4. Refine the MSA to improve quality
        input_msa = f"{msa_prefix}.fasta"
        refined_msa = os.path.join(output_dir, "opsin_refined.fasta")

        print("[INFO] Refining MSA with FoldMason...")
        run_refinemsa(structure_db, input_msa, refined_msa, refine_iters=1000)

        # 5. Generate LDDT data (avoid HTML format)
        # Instead of using msa2lddtreport, use msa2lddt to get LDDT scores directly
        from projects.opsin_analysis.foldmason_helpers import run_msa2lddt

        lddt_output = os.path.join(output_dir, "alignment_lddt_scores.txt")
        print("[INFO] Calculating LDDT scores...")
        try:
            # Run FoldMason's LDDT calculation and capture the output
            lddt_result = run_msa2lddt(structure_db, refined_msa)

            # Save the LDDT scores to a text file
            with open(lddt_output, 'w') as f:
                f.write(lddt_result)

            print("[INFO] LDDT scores saved to", lddt_output)
        except Exception as e:
            print(f"[WARNING] Error calculating LDDT scores: {str(e)}")
            lddt_output = None
            lddt_result = None

        # Process results as in original function...
        # Rest of the function implementation remains the same

        # Parse and analyze FoldMason results here...

        return {
            'processed_structures': processed_structures_complete,
            'foldmason_alignment': refined_msa,
            'foldmason_lddt_scores': lddt_output,
            'foldmason_lddt_result': lddt_result
        }

    except Exception as e:
        print(f"[ERROR] FoldMason alignment failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return data_dict
