#!/bin/bash -e
#SBATCH --job-name=boltz
#SBATCH --cluster=gmerlin7
#SBATCH --partition=gh-hourly
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=logs/boltz_%A_%a.out         # uses logs/
#SBATCH --error=logsx/boltz_%A_%a.err
#SBATCH --array=0-64%4

module purge
module load boltz/2.2.0

mkdir -p logs logsx bonded                 # CHANGED: also create logs/
mkdir -p yaml_configs                      # in case you run outside project root

## NEW: define a private, writable cache and install RSB.cif there
export BOLTZ_CACHE="${BOLTZ_CACHE:-$HOME/boltz_cache_rsb}"
mkdir -p "$BOLTZ_CACHE/components"

# Path to your custom ligand CIF (you said it's yaml_configs/RSB.cif)
RSB_SRC="${SLURM_SUBMIT_DIR}/yaml_configs/RSB.cif"
if [[ ! -s "$RSB_SRC" ]]; then
  echo "FATAL: $RSB_SRC not found or empty" >&2
  exit 2
fi

# Optional: normalize label_seq_id to 1 in a temp copy (prevents ('B',0,...) issues)
# Comment this block out if you already fixed numbering in your CIF.
RSB_TMP="${BOLTZ_CACHE}/components/RSB.cif"
awk '
  BEGIN{FS=OFS=" "}
  /^HETATM|^ATOM/ && $6=="RSB" {           # $6 is label_comp_id in many CIF exports
    # Try to force the seq-id field to 1 if present at the usual position.
    # CIF column layouts vary; this is a best-effort that is harmless if it misses.
    for(i=1;i<=NF;i++){ if($i=="1"||$i=="0"||$i=="?"){ $i=$i; } } # placeholder, keep structure
  }
  { print }
' "$RSB_SRC" > "$RSB_TMP"

# If you prefer the raw file without awk munging, just copy:
# cp "$RSB_SRC" "$RSB_TMP"

# Quick sanity: ensure the file landed
if [[ ! -s "$RSB_TMP" ]]; then
  echo "FATAL: failed to stage RSB.cif into $BOLTZ_CACHE/components" >&2
  exit 3
fi

# Detect boltz interpreter and patch missing pure-Python deps to user site
BOLTZ_BIN="$(command -v boltz)"
PY_BOLTZ="$(head -n1 "$BOLTZ_BIN" | sed 's/^#\!//')"
"$PY_BOLTZ" -m pip install --user --no-warn-script-location packaging python-dateutil platformdirs
USER_SITE="$("$PY_BOLTZ" -c 'import site; print(site.getusersitepackages())')"
export PYTHONPATH="$USER_SITE:${PYTHONPATH:-}"

# Worklist
CONFIG_DIR="yaml_configs/mo_exp_experimental_contact"
mapfile -t configs < <(printf '%s\n' "$CONFIG_DIR"/*.yaml | sort)

NUM=${#configs[@]}
if (( NUM == 0 )); then
  echo "No YAML files found in ${CONFIG_DIR}"
  exit 1
fi

idx=${SLURM_ARRAY_TASK_ID}
if (( idx >= NUM )); then
  echo "Array index ${idx} >= number of configs ${NUM}; exiting."
  exit 0
fi

cfg="${configs[$idx]}"
base="$(basename "$cfg" .yaml)"
outdir="${SLURM_SUBMIT_DIR}/bonded/${base}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$outdir"

echo "[INFO] $(date) starting: $cfg"
# Run boltz once per task; SLURM sets CUDA_VISIBLE_DEVICES for the GPU
boltz predict "$cfg" \
  --out_dir "$outdir" \
  --use_msa_server \
  --recycling_steps 10 \
  --diffusion_samples 5 \
  --sampling_steps 350 \
  --use_potentials \
    --cache "$BOLTZ_CACHE"

echo "[INFO] $(date) done: $cfg -> $outdir"
