#!/bin/bash
#SBATCH --job-name=GSA_quality
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/Grand_Saint_Antoine/00_scripts/02_quality.err
#SBATCH --output=/home/plstenge/Grand_Saint_Antoine/00_scripts/02_quality.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

# QC FastQC + MultiQC pour le projet Grand Saint Antoine uniquement
# Utilise une liste blanche stricte des dossiers du projet pour éviter
# d'inclure d'autres projets présents dans le même répertoire.

BASEDIR="/storage/groups/gdec/shared_paleo/default"
WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_qc_default"
THREADS="${THREADS:-16}"

mkdir -p "$WORKDIR"/{lists,fastqc,multiqc,logs}

GSA_DIRS=(
  1A-m1 2A-m1 3A-m1 4A-m1 6A-m1 NTC-A-m1
  1B-m1 2B-m1 3B-m1 4B-m1 6B-m1 NTC-B-m1
  1C-m1 2C-m1 3C-m1 4C-m1 6C-m1 NTC-C-m1
  Neg-m1
  1A-m2 2A-m2 3A-m2 4A-m2 6A-m2 NTC-A-m2
  1B-m2 2B-m2 3B-m2 4B-m2 6B-m2 NTC-B-m2
  1C-m2 2C-m2 3C-m2 4C-m2 6C-m2 NTC-C-m2
  Neg-m2
  1A-m3 2A-m3 3A-m3 4A-m3 6A-m3 NTC-A-m3
  1B-m3 2B-m3 3B-m3 4B-m3 6B-m3 NTC-B-m3
  1C-m3 2C-m3 3C-m3 4C-m3 6C-m3 NTC-C-m3
  Neg-m3
)

FASTQ_LIST="$WORKDIR/lists/gsa_fastq_files.txt"
: > "$FASTQ_LIST"

for sampledir in "${GSA_DIRS[@]}"; do
  d="$BASEDIR/$sampledir"
  if [[ ! -d "$d" ]]; then
    echo "ATTENTION: dossier absent: $d" | tee -a "$WORKDIR/logs/warnings.log" >&2
    continue
  fi

  find "$d" -maxdepth 1 -type f \( -iname "*_R1*.fastq.gz" -o -iname "*_R2*.fastq.gz" \) | sort >> "$FASTQ_LIST"
done

sort -u "$FASTQ_LIST" -o "$FASTQ_LIST"

NFILES=$(wc -l < "$FASTQ_LIST")
if [[ "$NFILES" -eq 0 ]]; then
  echo "ERREUR: aucun FASTQ trouvé pour le projet GSA." >&2
  exit 1
fi

echo "${NFILES} fichiers FASTQ détectés pour GSA"

if ! command -v fastqc >/dev/null 2>&1; then
  echo "ERREUR: fastqc non trouvé dans le PATH" >&2
  exit 1
fi

if ! command -v multiqc >/dev/null 2>&1; then
  echo "ERREUR: multiqc non trouvé dans le PATH" >&2
  exit 1
fi

# FastQC
fastqc \
  --threads "$THREADS" \
  --outdir "$WORKDIR/fastqc" \
  $(cat "$FASTQ_LIST") \
  2>&1 | tee "$WORKDIR/logs/fastqc.log"

# MultiQC
multiqc "$WORKDIR/fastqc" \
  --outdir "$WORKDIR/multiqc" \
  --filename multiqc_gsa_report \
  --force \
  2>&1 | tee "$WORKDIR/logs/multiqc.log"

echo
echo "QC terminé."
echo "FastQC  : $WORKDIR/fastqc"
echo "MultiQC : $WORKDIR/multiqc/multiqc_gsa_report.html"
echo "Liste des FASTQ utilisés : $FASTQ_LIST"
