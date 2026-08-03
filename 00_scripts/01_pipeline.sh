#!/bin/bash
#SBATCH --job-name=GSA
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline.err
#SBATCH --output=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

#############################################################################
# Pipeline QIIME2 - Analyse multi-marqueurs (pla, caf, V4V5)
# avec gestion des contrôles négatifs (NTC / Neg)
# Version corrigée complète
#############################################################################

########################
# 0. Paramètres généraux
########################
BASEDIR="/storage/groups/gdec/shared_paleo/E1739/filtermask"
WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_filtermask_pipeline"
THREADS=16
USE_ILLUMINA_ADAPTER_TRIM=0

mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata,logs}
cd "$WORKDIR"

# Primers spécifiques après les Ns
# m1: pla ; m2: caf1 ; m3: V4V5

declare -A PRIMER_F
declare -A PRIMER_R
PRIMER_F[m1]="GACTGGGTTCGGGCACATG"
PRIMER_R[m1]="AGACTTTGGCATTAGGTGTG"
PRIMER_F[m2]="aaccagcccgcatcactctta"
PRIMER_R[m2]="atcacccgcggcatctgta"
PRIMER_F[m3]="GTGYCAGCMGCCGCGGTAA"
PRIMER_R[m3]="CCGYCAATTYMTTTRAGTTT"

MARKERS=(m1 m2 m3)

# Adaptateurs Illumina universels
# Désactivés par défaut car cela avait généré un passage d'adaptateur vide dans q2-cutadapt.
# Active-les seulement si tu veux explicitement les enlever à cette étape.
ILLUMINA_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ILLUMINA_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

#####################################################
# 1. Construction des manifests (version précise GSA)
#    On ne prend QUE les dossiers explicitement attendus
#####################################################
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

for marker in "${MARKERS[@]}"; do
  manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
  printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$manifest"

  for sampledir in "${GSA_DIRS[@]}"; do
    [[ "$sampledir" == *-"$marker" ]] || continue
    d="$BASEDIR/$sampledir"

    if [[ ! -d "$d" ]]; then
      echo "ATTENTION: dossier absent: $d" >&2
      continue
    fi

    r1=$(find "$d" -maxdepth 1 -type f -iname "*_R1*.fastq.gz" | sort | head -n1)
    r2=$(find "$d" -maxdepth 1 -type f -iname "*_R2*.fastq.gz" | sort | head -n1)

    if [[ -z "$r1" || -z "$r2" ]]; then
      echo "ATTENTION: fastq manquants dans $d" >&2
      continue
    fi

    printf "%s\t%s\t%s\n" "$sampledir" "$r1" "$r2" >> "$manifest"
  done
done

#####################################################
# 2. Fichier de métadonnées (commun, avec statut contrôle)
#####################################################
meta="$WORKDIR/metadata/sample-metadata.tsv"
printf "sample-id\tsample\treplicate\tmarker\tsample-or-control\n" > "$meta"


for sampledir in "${GSA_DIRS[@]}"; do
  [[ "$sampledir" == *-"$marker" ]] || continue
  d="$BASEDIR/$sampledir"

  if [[ ! -d "$d" ]]; then
    echo "ATTENTION: dossier absent: $d" >&2
    continue
  fi

  r1=$(find "$d" -maxdepth 1 -type f -iname "*_R1*.fastq.gz" | sort | head -n1)
  r2=$(find "$d" -maxdepth 1 -type f -iname "*_R2*.fastq.gz" | sort | head -n1)

  if [[ -z "$r1" || -z "$r2" ]]; then
    echo "ATTENTION: fastq manquants dans $d" >&2
    continue
  fi

  printf "%s\t%s\t%s\n" "$sampledir" "$r1" "$r2" >> "$manifest"
done

echo "Metadata écrite : $meta"

#####################################################
# 3. Import QIIME2 + suppression primers/adaptateurs (cutadapt)
#####################################################
for marker in "${MARKERS[@]}"; do
  manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
  demux="$WORKDIR/imported/demux_${marker}.qza"
  trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
  trimmed_qzv="$WORKDIR/trimmed/trimmed_${marker}.qzv"
  cutadapt_stats="$WORKDIR/trimmed/cutadapt_${marker}_stats.qza"

  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$manifest" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$demux"

  if [[ "$USE_ILLUMINA_ADAPTER_TRIM" -eq 1 ]]; then
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences "$demux" \
      --p-cores "$THREADS" \
      --p-front-f "${PRIMER_F[$marker]}" \
      --p-front-r "${PRIMER_R[$marker]}" \
      --p-adapter-f "$ILLUMINA_FWD" \
      --p-adapter-r "$ILLUMINA_REV" \
      --p-match-read-wildcards \
      --p-discard-untrimmed \
      --p-no-indels \
      --o-trimmed-sequences "$trimmed" \
      --o-stats "$cutadapt_stats" \
      --verbose \
      2>&1 | tee "$WORKDIR/logs/cutadapt_${marker}.log"
  else
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences "$demux" \
      --p-cores "$THREADS" \
      --p-front-f "${PRIMER_F[$marker]}" \
      --p-front-r "${PRIMER_R[$marker]}" \
      --p-match-read-wildcards \
      --p-discard-untrimmed \
      --p-no-indels \
      --o-trimmed-sequences "$trimmed" \
      --o-stats "$cutadapt_stats" \
      --verbose \
      2>&1 | tee "$WORKDIR/logs/cutadapt_${marker}.log"
  fi

  [[ -f "$trimmed" ]] || { echo "ERREUR: cutadapt n'a pas créé $trimmed" >&2; exit 1; }

  qiime demux summarize \
    --i-data "$trimmed" \
    --o-visualization "$trimmed_qzv"
done

#####################################################
# 4. DADA2 (denoise-paired) par marqueur
#    IMPORTANT : inspecter les .qzv de l'étape 3 pour ajuster si besoin
#####################################################
declare -A TRUNC_F
declare -A TRUNC_R
TRUNC_F[m1]=0
TRUNC_R[m1]=0
TRUNC_F[m2]=0
TRUNC_R[m2]=0
TRUNC_F[m3]=280
TRUNC_R[m3]=220

for marker in "${MARKERS[@]}"; do
  trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
  table="$WORKDIR/dada2/table_${marker}.qza"
  repseqs="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  stats="$WORKDIR/dada2/stats_${marker}.qza"

  [[ -f "$trimmed" ]] || { echo "ERREUR: fichier absent $trimmed" >&2; exit 1; }

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$trimmed" \
    --p-trunc-len-f "${TRUNC_F[$marker]}" \
    --p-trunc-len-r "${TRUNC_R[$marker]}" \
    --p-n-threads "$THREADS" \
    --o-table "$table" \
    --o-representative-sequences "$repseqs" \
    --o-denoising-stats "$stats" \
    --verbose

  qiime feature-table summarize \
    --i-table "$table" \
    --m-sample-metadata-file "$meta" \
    --o-visualization "$WORKDIR/dada2/table_${marker}.qzv"

  qiime feature-table tabulate-seqs \
    --i-data "$repseqs" \
    --o-visualization "$WORKDIR/dada2/rep-seqs_${marker}.qzv"
done

#####################################################
# 5. Taxonomie
#    - V4V5 (m3) : classificateur SILVA téléchargé si absent
#    - pla / caf (m1, m2) : base NCBI construite si absente
#####################################################
SILVA_URL="https://data.qiime2.org/2024.5/common/silva-138-99-nb-classifier.qza"
CLASSIFIER="$WORKDIR/taxonomy/silva-138-99-nb-classifier.qza"

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "Classificateur SILVA absent -> téléchargement en cours..."
  wget -c -O "$CLASSIFIER" "$SILVA_URL"
  if ! unzip -l "$CLASSIFIER" >/dev/null 2>&1; then
    echo "ERREUR : le téléchargement du classificateur SILVA a échoué (fichier invalide)." >&2
    rm -f "$CLASSIFIER"
  fi
fi

if [[ -f "$CLASSIFIER" ]]; then
  qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$WORKDIR/dada2/rep-seqs_m3.qza" \
    --p-n-jobs "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m3.qza"
else
  echo "ATTENTION : classification V4V5 impossible, classificateur SILVA introuvable." >&2
fi

declare -A NCBI_QUERY
NCBI_QUERY[m1]='(Yersinia pestis[Organism] AND pla[Gene]) AND (plasminogen activator[Title] OR pla[Title])'
NCBI_QUERY[m2]='(Yersinia pestis[Organism] OR Enterobacteriaceae[Organism]) AND (caf1[Gene] OR "capsular antigen F1"[Title] OR caf1[Title])'

for marker in m1 m2; do
  REF_SEQS="$WORKDIR/taxonomy/ref-seqs_${marker}.qza"
  REF_TAX="$WORKDIR/taxonomy/ref-taxonomy_${marker}.qza"

  if [[ ! -f "$REF_SEQS" || ! -f "$REF_TAX" ]]; then
    echo "Base de référence ${marker} absente -> construction automatique depuis NCBI..."
    RAW_SEQS="$WORKDIR/taxonomy/ncbi-${marker}-seqs-raw.qza"
    RAW_TAX="$WORKDIR/taxonomy/ncbi-${marker}-tax-raw.qza"
    CULLED_SEQS="$WORKDIR/taxonomy/ncbi-${marker}-seqs-culled.qza"
    DISCARDED_SEQS="$WORKDIR/taxonomy/ncbi-${marker}-seqs-discarded.qza"

    qiime rescript get-ncbi-data \
      --p-query "${NCBI_QUERY[$marker]}" \
      --p-n-jobs 3 \
      --o-sequences "$RAW_SEQS" \
      --o-taxonomy "$RAW_TAX" \
      --verbose

    qiime rescript cull-seqs \
      --i-sequences "$RAW_SEQS" \
      --p-num-degenerates 5 \
      --p-homopolymer-length 12 \
      --o-clean-sequences "$CULLED_SEQS"

    qiime rescript filter-seqs-length \
      --i-sequences "$CULLED_SEQS" \
      --p-global-min 80 \
      --p-global-max 2000 \
      --o-filtered-seqs "$REF_SEQS" \
      --o-discarded-seqs "$DISCARDED_SEQS"

    cp "$RAW_TAX" "$REF_TAX"
  fi

  if [[ -f "$REF_SEQS" && -f "$REF_TAX" ]]; then
    qiime feature-classifier classify-consensus-vsearch \
      --i-query "$WORKDIR/dada2/rep-seqs_${marker}.qza" \
      --i-reference-reads "$REF_SEQS" \
      --i-reference-taxonomy "$REF_TAX" \
      --p-threads "$THREADS" \
      --o-classification "$WORKDIR/taxonomy/taxonomy_${marker}.qza" \
      --o-search-results "$WORKDIR/taxonomy/search_${marker}.qza"
  else
    echo "ATTENTION : classification ${marker} impossible (base de référence indisponible)." >&2
  fi
done

#####################################################
# 6. Inspection des contrôles AVANT nettoyage
#####################################################
for marker in "${MARKERS[@]}"; do
  taxo="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
  [ -f "$taxo" ] || continue

  qiime feature-table filter-samples \
    --i-table "$WORKDIR/dada2/table_${marker}.qza" \
    --m-metadata-file "$meta" \
    --p-where "[marker]='${marker}' AND [sample-or-control]='control'" \
    --o-filtered-table "$WORKDIR/decontam/controls_table_${marker}.qza"

  qiime taxa barplot \
    --i-table "$WORKDIR/decontam/controls_table_${marker}.qza" \
    --i-taxonomy "$taxo" \
    --m-metadata-file "$meta" \
    --o-visualization "$WORKDIR/decontam/controls_barplot_${marker}.qzv"
done

#####################################################
# 7. Décontamination (decontam, méthode prevalence)
#####################################################
for marker in "${MARKERS[@]}"; do
  table="$WORKDIR/dada2/table_${marker}.qza"
  repseqs="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  [ -f "$table" ] || continue

  qiime quality-control decontam-identify \
    --i-table "$table" \
    --m-metadata-file "$meta" \
    --p-method prevalence \
    --p-prev-control-column sample-or-control \
    --p-prev-control-indicator control \
    --o-decontam-scores "$WORKDIR/decontam/decontam-scores_${marker}.qza"

  qiime quality-control decontam-score-viz \
    --i-decontam-scores "$WORKDIR/decontam/decontam-scores_${marker}.qza" \
    --i-table "$table" \
    --p-threshold 0.1 \
    --o-visualization "$WORKDIR/decontam/decontam-scoreviz_${marker}.qzv"

  qiime quality-control decontam-remove \
    --i-decontam-scores "$WORKDIR/decontam/decontam-scores_${marker}.qza" \
    --i-table "$table" \
    --i-rep-seqs "$repseqs" \
    --p-threshold 0.1 \
    --o-filtered-table "$WORKDIR/decontam/table_${marker}_decontam.qza" \
    --o-filtered-rep-seqs "$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"
done

#####################################################
# 8. Retrait des échantillons contrôles eux-mêmes de la table finale
#####################################################
for marker in "${MARKERS[@]}"; do
  table_dc="$WORKDIR/decontam/table_${marker}_decontam.qza"
  repseqs_dc="$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"
  [ -f "$table_dc" ] || continue

  qiime feature-table filter-samples \
    --i-table "$table_dc" \
    --m-metadata-file "$meta" \
    --p-where "[marker]='${marker}' AND [sample-or-control]='sample'" \
    --o-filtered-table "$WORKDIR/decontam/table_${marker}_final.qza"

  qiime feature-table filter-seqs \
    --i-data "$repseqs_dc" \
    --i-table "$WORKDIR/decontam/table_${marker}_final.qza" \
    --o-filtered-data "$WORKDIR/decontam/rep-seqs_${marker}_final.qza"
done

#####################################################
# 9. Table taxonomique finale + exports (pas de raréfaction)
#####################################################
for marker in "${MARKERS[@]}"; do
  table_final="$WORKDIR/decontam/table_${marker}_final.qza"
  taxo="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
  [ -f "$table_final" ] && [ -f "$taxo" ] || continue

  qiime taxa barplot \
    --i-table "$table_final" \
    --i-taxonomy "$taxo" \
    --m-metadata-file "$meta" \
    --o-visualization "$WORKDIR/exports/taxa-barplot_${marker}_final.qzv"

  mkdir -p "$WORKDIR/exports/${marker}_final"
  qiime tools export --input-path "$table_final" --output-path "$WORKDIR/exports/${marker}_final"
  qiime tools export --input-path "$taxo" --output-path "$WORKDIR/exports/${marker}_final"

  biom add-metadata \
    -i "$WORKDIR/exports/${marker}_final/feature-table.biom" \
    -o "$WORKDIR/exports/${marker}_final/feature-table-with-tax.biom" \
    --observation-metadata-fp "$WORKDIR/exports/${marker}_final/taxonomy.tsv" \
    --observation-header "OTUID,taxonomy,confidence" \
    --sc-separated taxonomy

  biom convert \
    -i "$WORKDIR/exports/${marker}_final/feature-table-with-tax.biom" \
    -o "$WORKDIR/exports/${marker}_final/ASV_table_${marker}_taxonomy.tsv" \
    --to-tsv --header-key taxonomy

done

echo "Pipeline terminé."
