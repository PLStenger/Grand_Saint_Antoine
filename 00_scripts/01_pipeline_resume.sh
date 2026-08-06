#!/bin/bash
#SBATCH --job-name=GSA
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline_resume.err
#SBATCH --output=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline_resume.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

#############################################################################
# Pipeline QIIME2 - reprise / finalisation
# GSA - pla (m1), caf1 (m2), V4V5 (m3)
#############################################################################

########################
# 0. Paramètres généraux
########################
BASEDIR="/storage/groups/gdec/shared_paleo/E1739/filtermask"
WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_filtermask_pipeline"
THREADS=16

mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata,logs,tmp}
cd "$WORKDIR"

META="$WORKDIR/metadata/sample-metadata.tsv"

########################
# 0bis. Vérifs minimales
########################
echo "=== Vérification des fichiers existants ==="

for f in \
  "$WORKDIR/trimmed/trimmed_m1.qza" \
  "$WORKDIR/trimmed/trimmed_m2.qza" \
  "$WORKDIR/trimmed/trimmed_m3.qza" \
  "$WORKDIR/taxonomy/ref-seqs_m1.qza" \
  "$WORKDIR/taxonomy/ref-taxonomy_m1.qza" \
  "$WORKDIR/taxonomy/ref-seqs_m2.qza" \
  "$WORKDIR/taxonomy/ref-taxonomy_m2.qza"
do
  [[ -f "$f" ]] || { echo "ERREUR: fichier attendu absent: $f" >&2; exit 1; }
done

[[ -f "$META" ]] || { echo "ERREUR: metadata absente: $META" >&2; exit 1; }

########################
# 1. Étapes déjà faites : on ne relance pas
########################
# On garde la trace de ce qui a déjà été exécuté, mais on ne relance pas.
#
# # Construction des manifests
# # Import QIIME2
# # Cutadapt sur m1/m2/m3
# # Construction des bases NCBI pour m1/m2
# #
# # Fichiers déjà présents observés :
# # trimmed/trimmed_m1.qza
# # trimmed/trimmed_m2.qza
# # trimmed/trimmed_m3.qza
# # taxonomy/ref-seqs_m1.qza
# # taxonomy/ref-taxonomy_m1.qza
# # taxonomy/ref-seqs_m2.qza
# # taxonomy/ref-taxonomy_m2.qza
# # taxonomy/taxonomy_m1.qza
# # taxonomy/search_m1.qza
# # taxonomy/search_m2.qza
# #
# # Donc on reprend à partir du débruitage / classification / export.

########################
# 2. Paramètres DADA2
########################

# m1 : déjà OK, on ne relance pas
# m2 : rerun paired-end avec paramètres relâchés
# m3 : rerun en single-end sur R1 seulement

M2_TRUNC_F=0
M2_TRUNC_R=0
M2_MAXEE_F=8.0
M2_MAXEE_R=8.0
M2_TRIMLEFT_F=0
M2_TRIMLEFT_R=0

# m3 : single-end sur forward uniquement
M3_TRUNC_LEN=220
M3_MAXEE=8.0
M3_TRIMLEFT=0

########################
# 3. m1 - résumés si besoin
########################
echo "=== m1 : vérification des artefacts existants ==="

if [[ -f "$WORKDIR/dada2/table_m1.qza" && -f "$WORKDIR/dada2/rep-seqs_m1.qza" ]]; then
  echo "m1 DADA2 déjà présent."
else
  echo "ATTENTION: m1 n'a pas de table/rep-seqs, ce n'était pas attendu." >&2
fi

if [[ -f "$WORKDIR/dada2/table_m1.qza" ]]; then
  qiime feature-table summarize \
    --i-table "$WORKDIR/dada2/table_m1.qza" \
    --m-metadata-file "$META" \
    --o-feature-frequencies "$WORKDIR/dada2/feature-freq_m1.qza" \
    --o-sample-frequencies "$WORKDIR/dada2/sample-freq_m1.qza" \
    --o-summary "$WORKDIR/dada2/table_m1.qzv" \
    || true
fi

if [[ -f "$WORKDIR/dada2/rep-seqs_m1.qza" ]]; then
  qiime feature-table tabulate-seqs \
    --i-data "$WORKDIR/dada2/rep-seqs_m1.qza" \
    --o-visualization "$WORKDIR/dada2/rep-seqs_m1.qzv" \
    || true
fi

########################
# 4. m2 - rerun DADA2 paired-end
########################
echo "=== m2 : rerun DADA2 paired-end ==="

# On sauvegarde les anciens résultats s'ils existent
for f in \
  "$WORKDIR/dada2/table_m2.qza" \
  "$WORKDIR/dada2/rep-seqs_m2.qza" \
  "$WORKDIR/dada2/stats_m2.qza" \
  "$WORKDIR/dada2/base-transition-stats_m2.qza" \
  "$WORKDIR/dada2/table_m2.qzv" \
  "$WORKDIR/dada2/rep-seqs_m2.qzv"
do
  [[ -f "$f" ]] && mv "$f" "${f}.bak_$(date +%Y%m%d_%H%M%S)" || true
done

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$WORKDIR/trimmed/trimmed_m2.qza" \
  --p-trunc-len-f "$M2_TRUNC_F" \
  --p-trunc-len-r "$M2_TRUNC_R" \
  --p-trim-left-f "$M2_TRIMLEFT_F" \
  --p-trim-left-r "$M2_TRIMLEFT_R" \
  --p-max-ee-f "$M2_MAXEE_F" \
  --p-max-ee-r "$M2_MAXEE_R" \
  --p-n-threads "$THREADS" \
  --p-n-reads-learn 100000 \
  --o-table "$WORKDIR/dada2/table_m2.qza" \
  --o-representative-sequences "$WORKDIR/dada2/rep-seqs_m2.qza" \
  --o-denoising-stats "$WORKDIR/dada2/stats_m2.qza" \
  --o-base-transition-stats "$WORKDIR/dada2/base-transition-stats_m2.qza" \
  --verbose \
  2>&1 | tee "$WORKDIR/logs/dada2_m2_rerun.log"

if [[ -f "$WORKDIR/dada2/table_m2.qza" ]]; then
  qiime feature-table summarize \
    --i-table "$WORKDIR/dada2/table_m2.qza" \
    --m-metadata-file "$META" \
    --o-feature-frequencies "$WORKDIR/dada2/feature-freq_m2.qza" \
    --o-sample-frequencies "$WORKDIR/dada2/sample-freq_m2.qza" \
    --o-summary "$WORKDIR/dada2/table_m2.qzv"

  qiime feature-table tabulate-seqs \
    --i-data "$WORKDIR/dada2/rep-seqs_m2.qza" \
    --o-visualization "$WORKDIR/dada2/rep-seqs_m2.qzv"
else
  echo "ERREUR: m2 DADA2 n'a pas produit de table." >&2
fi

########################
# 5. m3 - extraction forward et DADA2 single-end
########################
echo "=== m3 : extraction forward puis DADA2 single-end ==="

# On crée un manifest single-end à partir du manifest m3 paired existant
M3_MANIFEST_PAIRED="$WORKDIR/manifests/manifest_m3.tsv"
M3_MANIFEST_FWD="$WORKDIR/manifests/manifest_m3_forward.tsv"

[[ -f "$M3_MANIFEST_PAIRED" ]] || { echo "ERREUR: manifest m3 absent." >&2; exit 1; }

printf "sample-id\tabsolute-filepath\n" > "$M3_MANIFEST_FWD"
tail -n +2 "$M3_MANIFEST_PAIRED" | awk -F '\t' '{print $1 "\t" $2}' >> "$M3_MANIFEST_FWD"

# Sauvegarde si déjà présent
[[ -f "$WORKDIR/imported/demux_m3_forward.qza" ]] && mv "$WORKDIR/imported/demux_m3_forward.qza" "$WORKDIR/imported/demux_m3_forward.qza.bak_$(date +%Y%m%d_%H%M%S)" || true

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$M3_MANIFEST_FWD" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path "$WORKDIR/imported/demux_m3_forward.qza"

# Sauvegarde anciens résultats m3 si présents
for f in \
  "$WORKDIR/dada2/table_m3.qza" \
  "$WORKDIR/dada2/rep-seqs_m3.qza" \
  "$WORKDIR/dada2/stats_m3.qza" \
  "$WORKDIR/dada2/base-transition-stats_m3.qza" \
  "$WORKDIR/dada2/table_m3.qzv" \
  "$WORKDIR/dada2/rep-seqs_m3.qzv"
do
  [[ -f "$f" ]] && mv "$f" "${f}.bak_$(date +%Y%m%d_%H%M%S)" || true
done

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "$WORKDIR/imported/demux_m3_forward.qza" \
  --p-trunc-len "$M3_TRUNC_LEN" \
  --p-trim-left "$M3_TRIMLEFT" \
  --p-max-ee "$M3_MAXEE" \
  --p-n-threads "$THREADS" \
  --p-n-reads-learn 100000 \
  --o-table "$WORKDIR/dada2/table_m3.qza" \
  --o-representative-sequences "$WORKDIR/dada2/rep-seqs_m3.qza" \
  --o-denoising-stats "$WORKDIR/dada2/stats_m3.qza" \
  --o-base-transition-stats "$WORKDIR/dada2/base-transition-stats_m3.qza" \
  --verbose \
  2>&1 | tee "$WORKDIR/logs/dada2_m3_single.log"

if [[ -f "$WORKDIR/dada2/table_m3.qza" ]]; then
  qiime feature-table summarize \
    --i-table "$WORKDIR/dada2/table_m3.qza" \
    --m-metadata-file "$META" \
    --o-feature-frequencies "$WORKDIR/dada2/feature-freq_m3.qza" \
    --o-sample-frequencies "$WORKDIR/dada2/sample-freq_m3.qza" \
    --o-summary "$WORKDIR/dada2/table_m3.qzv"

  qiime feature-table tabulate-seqs \
    --i-data "$WORKDIR/dada2/rep-seqs_m3.qza" \
    --o-visualization "$WORKDIR/dada2/rep-seqs_m3.qzv"
else
  echo "ERREUR: m3 DADA2 single n'a pas produit de table." >&2
fi

########################
# 6. Taxonomie
########################
echo "=== Taxonomie ==="

# m1 et m2 : bases déjà faites ; relancer classification si rep-seqs présents
for marker in m1 m2; do
  REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  REF_SEQS="$WORKDIR/taxonomy/ref-seqs_${marker}.qza"
  REF_TAX="$WORKDIR/taxonomy/ref-taxonomy_${marker}.qza"
  OUT_TAX="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
  OUT_SEARCH="$WORKDIR/taxonomy/search_${marker}.qza"

  if [[ -f "$REP" && -f "$REF_SEQS" && -f "$REF_TAX" ]]; then
    qiime feature-classifier classify-consensus-vsearch \
      --i-query "$REP" \
      --i-reference-reads "$REF_SEQS" \
      --i-reference-taxonomy "$REF_TAX" \
      --p-threads "$THREADS" \
      --o-classification "$OUT_TAX" \
      --o-search-results "$OUT_SEARCH"
  else
    echo "ATTENTION: taxonomie ${marker} non relancée (fichiers manquants)." >&2
  fi
done

# m3 : nécessite un classificateur SILVA 515F-926R déjà entraîné / disponible
M3_CLASSIFIER="$WORKDIR/taxonomy/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [[ -f "$M3_CLASSIFIER" && -f "$WORKDIR/dada2/rep-seqs_m3.qza" ]]; then
  qiime feature-classifier classify-sklearn \
    --i-classifier "$M3_CLASSIFIER" \
    --i-reads "$WORKDIR/dada2/rep-seqs_m3.qza" \
    --p-n-jobs "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m3.qza"
else
  echo "ATTENTION : taxonomie m3 non faite (classificateur SILVA 515F-926R absent ou rep-seqs_m3 absent)." >&2
fi

########################
# 7. Contrôles négatifs : inspection visuelle
########################
echo "=== Inspection des contrôles ==="

for marker in m1 m2 m3; do
  TABLE="$WORKDIR/dada2/table_${marker}.qza"
  TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
  CTRL_TABLE="$WORKDIR/decontam/controls_table_${marker}.qza"
  CTRL_BAR="$WORKDIR/decontam/controls_barplot_${marker}.qzv"

  [[ -f "$TABLE" ]] || continue
  [[ -f "$TAXO" ]] || continue

  # On filtre directement les contrôles via metadata
  qiime feature-table filter-samples \
    --i-table "$TABLE" \
    --m-metadata-file "$META" \
    --p-where "[sample-or-control]='control' AND [marker]='${marker}'" \
    --o-filtered-table "$CTRL_TABLE" \
    || true

  if [[ -f "$CTRL_TABLE" ]]; then
    qiime taxa barplot \
      --i-table "$CTRL_TABLE" \
      --i-taxonomy "$TAXO" \
      --m-metadata-file "$META" \
      --o-visualization "$CTRL_BAR" \
      || true
  fi
done

########################
# 8. Decontam : score seulement
########################
echo "=== Decontam : identification ==="

for marker in m1 m2 m3; do
  TABLE="$WORKDIR/dada2/table_${marker}.qza"
  REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  SCORE="$WORKDIR/decontam/decontam-scores_${marker}.qza"
  SCORE_VIZ="$WORKDIR/decontam/decontam-scoreviz_${marker}.qzv"

  [[ -f "$TABLE" ]] || continue

  qiime quality-control decontam-identify \
    --i-table "$TABLE" \
    --m-metadata-file "$META" \
    --p-method prevalence \
    --p-prev-control-column sample-or-control \
    --p-prev-control-indicator control \
    --o-decontam-scores "$SCORE" \
    || true

  if [[ -f "$SCORE" ]]; then
    if [[ -f "$REP" ]]; then
      qiime quality-control decontam-score-viz \
        --i-decontam-scores "$SCORE" \
        --i-table "$TABLE" \
        --i-rep-seqs "$REP" \
        --p-threshold 0.1 \
        --o-visualization "$SCORE_VIZ" \
        || true
    else
      qiime quality-control decontam-score-viz \
        --i-decontam-scores "$SCORE" \
        --i-table "$TABLE" \
        --p-threshold 0.1 \
        --o-visualization "$SCORE_VIZ" \
        || true
    fi
  fi
done

########################
# 9. Retrait des contaminants sans decontam-remove
########################
echo "=== Filtrage des contaminants ==="

for marker in m1 m2 m3; do
  TABLE="$WORKDIR/dada2/table_${marker}.qza"
  REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  SCORE="$WORKDIR/decontam/decontam-scores_${marker}.qza"

  TABLE_DC="$WORKDIR/decontam/table_${marker}_decontam.qza"
  REP_DC="$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"

  TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
  REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"

  [[ -f "$TABLE" ]] || continue
  [[ -f "$REP" ]] || continue
  [[ -f "$SCORE" ]] || continue

  qiime feature-table filter-features \
    --i-table "$TABLE" \
    --m-metadata-file "$SCORE" \
    --p-where '[p] > 0.1 OR [p] IS NULL' \
    --o-filtered-table "$TABLE_DC" \
    || true

  if [[ -f "$TABLE_DC" ]]; then
    qiime feature-table filter-seqs \
      --i-data "$REP" \
      --i-table "$TABLE_DC" \
      --o-filtered-data "$REP_DC" \
      || true
  fi

  # retrait des échantillons contrôles eux-mêmes dans la table finale
  if [[ -f "$TABLE_DC" ]]; then
    qiime feature-table filter-samples \
      --i-table "$TABLE_DC" \
      --m-metadata-file "$META" \
      --p-where "[sample-or-control]='sample' AND [marker]='${marker}'" \
      --o-filtered-table "$TABLE_FINAL" \
      || true
  fi

  if [[ -f "$TABLE_FINAL" && -s "$TABLE_FINAL" ]]; then
    qiime feature-table filter-seqs \
      --i-data "$REP_DC" \
      --i-table "$TABLE_FINAL" \
      --o-filtered-data "$REP_FINAL" \
      || true
  fi
done

########################
# 10. Exports finaux
########################
echo "=== Exports finaux ==="

for marker in m1 m2 m3; do
  TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"

  TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
  REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"

  TABLE_RAW="$WORKDIR/dada2/table_${marker}.qza"
  REP_RAW="$WORKDIR/dada2/rep-seqs_${marker}.qza"

  EXPORT_DIR="$WORKDIR/exports/${marker}_final"

  mkdir -p "$EXPORT_DIR"

  # Si la table finale existe, on l'utilise ; sinon fallback sur la table brute
  TABLE_TO_USE="$TABLE_RAW"
  REP_TO_USE="$REP_RAW"

  if [[ -f "$TABLE_FINAL" && -f "$REP_FINAL" ]]; then
    TABLE_TO_USE="$TABLE_FINAL"
    REP_TO_USE="$REP_FINAL"
  fi

  [[ -f "$TABLE_TO_USE" ]] || continue
  [[ -f "$REP_TO_USE" ]] || continue

  qiime feature-table summarize \
    --i-table "$TABLE_TO_USE" \
    --m-metadata-file "$META" \
    --o-feature-frequencies "$EXPORT_DIR/feature-frequencies_${marker}.qza" \
    --o-sample-frequencies "$EXPORT_DIR/sample-frequencies_${marker}.qza" \
    --o-summary "$EXPORT_DIR/table_${marker}.qzv" \
    || true

  qiime feature-table tabulate-seqs \
    --i-data "$REP_TO_USE" \
    --o-visualization "$EXPORT_DIR/rep-seqs_${marker}.qzv" \
    || true

  if [[ -f "$TAXO" ]]; then
    qiime taxa barplot \
      --i-table "$TABLE_TO_USE" \
      --i-taxonomy "$TAXO" \
      --m-metadata-file "$META" \
      --o-visualization "$EXPORT_DIR/taxa-barplot_${marker}.qzv" \
      || true
  fi

  qiime tools export \
    --input-path "$TABLE_TO_USE" \
    --output-path "$EXPORT_DIR/table_export"

  qiime tools export \
    --input-path "$REP_TO_USE" \
    --output-path "$EXPORT_DIR/repseq_export"

  if [[ -f "$TAXO" ]]; then
    qiime tools export \
      --input-path "$TAXO" \
      --output-path "$EXPORT_DIR/taxonomy_export"
  fi

  if [[ -f "$EXPORT_DIR/table_export/feature-table.biom" && -f "$EXPORT_DIR/taxonomy_export/taxonomy.tsv" ]]; then
    biom add-metadata \
      -i "$EXPORT_DIR/table_export/feature-table.biom" \
      -o "$EXPORT_DIR/feature-table-with-tax.biom" \
      --observation-metadata-fp "$EXPORT_DIR/taxonomy_export/taxonomy.tsv" \
      --observation-header "OTUID,taxonomy,confidence" \
      --sc-separated taxonomy

    biom convert \
      -i "$EXPORT_DIR/feature-table-with-tax.biom" \
      -o "$EXPORT_DIR/ASV_table_${marker}_taxonomy.tsv" \
      --to-tsv --header-key taxonomy
  fi
done

echo "=== Pipeline de reprise terminé ==="
