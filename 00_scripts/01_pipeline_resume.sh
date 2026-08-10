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
# Pipeline QIIME2 - reprise v3 / finalisation
# GSA - pla (m1), caf1 (m2), V4V5 (m3)
#
# Corrections apportees suite au run precedent (01_pipeline_resume.out) :
#   1) m2 DADA2 : "No features remain after denoising" meme avec maxEE=8.
#      -> passage en pooling pseudo + chimera_method none pour recuperer
#         des ASV, puis on retente avec chimera consensus en secours.
#   2) m3 taxonomie : le classificateur SILVA est en realite dans
#      /home/plstenge/Grand_Saint_Antoine/taxonomy/ et non dans
#      02_filtermask_pipeline/taxonomy/ -> chemin corrige.
#   3) m1 100% Unassigned : classify-consensus-vsearch avec les seuils
#      par defaut (perc-identity 0.97, min-consensus 0.51) est trop strict
#      pour la base pla reduite -> seuils relaches + diagnostic sur
#      search_m1.qza (nombre de hits bruts) avant de conclure.
#############################################################################

########################
# 0. Parametres generaux
########################
BASEDIR="/storage/groups/gdec/shared_paleo/E1739/filtermask"
WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_filtermask_pipeline"
TAXODIR_GLOBAL="/home/plstenge/Grand_Saint_Antoine/taxonomy"
THREADS=16

mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata,logs,tmp}
cd "$WORKDIR"

META="$WORKDIR/metadata/sample-metadata.tsv"

########################
# 0bis. Verifs minimales
########################
echo "=== Verification des fichiers existants ==="

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
# 1. Etapes deja faites : on NE relance PAS
########################
# Tout ce bloc est laisse en commentaire pour trace, ne pas decommenter.
#
# # Construction des manifests (m1/m2/m3)
# # Import QIIME2 (m1/m2/m3)
# # Cutadapt sur m1/m2/m3            -> trimmed/trimmed_m{1,2,3}.qza
# # Construction des bases NCBI m1/m2 -> taxonomy/ref-seqs_m{1,2}.qza,
# #                                      taxonomy/ref-taxonomy_m{1,2}.qza
# # DADA2 m1 (deja OK, ne pas retoucher) -> dada2/table_m1.qza, rep-seqs_m1.qza
# # DADA2 m3 single-end (deja OK)        -> dada2/table_m3.qza, rep-seqs_m3.qza
# # Taxonomie m1 (classify-consensus-vsearch, deja fait mais tout Unassigned,
# #               on va relancer avec seuils relaches, cf etape 6)
#
# Donc dans ce script v3 on ne relance QUE :
#   - m2 DADA2 (avec strategie alternative)
#   - m3 taxonomie (chemin classificateur corrige)
#   - m1 taxonomie (seuils relaches + diagnostic)
#   - decontam + exports pour tout ce qui est disponible

########################
# 2. m1 - deja present, on ne relance pas DADA2
########################
echo "=== m1 : DADA2 deja present, aucune action ==="
if [[ ! -f "$WORKDIR/dada2/table_m1.qza" || ! -f "$WORKDIR/dada2/rep-seqs_m1.qza" ]]; then
  echo "ERREUR: m1 devrait deja avoir table/rep-seqs, verifier l'etat du dossier." >&2
fi

########################
# 3. m2 - rerun DADA2 avec strategie alternative
########################
echo "=== m2 : rerun DADA2 paired-end (pooling pseudo, sans puis avec chimere) ==="

for f in \
  "$WORKDIR/dada2/table_m2.qza" \
  "$WORKDIR/dada2/rep-seqs_m2.qza" \
  "$WORKDIR/dada2/stats_m2.qza" \
  "$WORKDIR/dada2/base-transition-stats_m2.qza" \
  "$WORKDIR/dada2/table_m2.qzv" \
  "$WORKDIR/dada2/rep-seqs_m2.qzv"
do
  [[ -f "$f" ]] && mv "$f" "${f}.bak_$(date +%Y%m%d_%H%M%S)"
done

M2_TRUNC_F=0
M2_TRUNC_R=0
M2_MAXEE_F=10.0
M2_MAXEE_R=10.0
M2_TRIMLEFT_F=0
M2_TRIMLEFT_R=0

# Tentative 1 : pooling pseudo, sans detection de chimeres
# (le pooling "independent" par echantillon est trop strict quand
#  chaque echantillon a tres peu de reads ; "pseudo" partage l'info
#  entre echantillons pour stabiliser l'apprentissage des erreurs)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$WORKDIR/trimmed/trimmed_m2.qza" \
  --p-trunc-len-f "$M2_TRUNC_F" \
  --p-trunc-len-r "$M2_TRUNC_R" \
  --p-trim-left-f "$M2_TRIMLEFT_F" \
  --p-trim-left-r "$M2_TRIMLEFT_R" \
  --p-max-ee-f "$M2_MAXEE_F" \
  --p-max-ee-r "$M2_MAXEE_R" \
  --p-pooling-method pseudo \
  --p-chimera-method none \
  --p-n-threads "$THREADS" \
  --p-n-reads-learn 100000 \
  --o-table "$WORKDIR/dada2/table_m2.qza" \
  --o-representative-sequences "$WORKDIR/dada2/rep-seqs_m2.qza" \
  --o-denoising-stats "$WORKDIR/dada2/stats_m2.qza" \
  --o-base-transition-stats "$WORKDIR/dada2/base-transition-stats_m2.qza" \
  --verbose \
  2>&1 | tee "$WORKDIR/logs/dada2_m2_pseudo_nochimera.log"

if [[ ! -f "$WORKDIR/dada2/table_m2.qza" ]]; then
  echo "Tentative 1 (pseudo, sans chimere) a echoue, essai 2 : pooling pseudo + chimere consensus" >&2

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$WORKDIR/trimmed/trimmed_m2.qza" \
    --p-trunc-len-f "$M2_TRUNC_F" \
    --p-trunc-len-r "$M2_TRUNC_R" \
    --p-trim-left-f "$M2_TRIMLEFT_F" \
    --p-trim-left-r "$M2_TRIMLEFT_R" \
    --p-max-ee-f "$M2_MAXEE_F" \
    --p-max-ee-r "$M2_MAXEE_R" \
    --p-pooling-method pseudo \
    --p-chimera-method consensus \
    --p-n-threads "$THREADS" \
    --p-n-reads-learn 100000 \
    --o-table "$WORKDIR/dada2/table_m2.qza" \
    --o-representative-sequences "$WORKDIR/dada2/rep-seqs_m2.qza" \
    --o-denoising-stats "$WORKDIR/dada2/stats_m2.qza" \
    --o-base-transition-stats "$WORKDIR/dada2/base-transition-stats_m2.qza" \
    --verbose \
    2>&1 | tee "$WORKDIR/logs/dada2_m2_pseudo_consensus.log"
fi

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

  echo "m2 : DADA2 a produit une table (voir logs pour savoir si sans ou avec detection de chimeres)."
else
  echo "ERREUR: m2 DADA2 n'a produit aucune table meme avec pooling pseudo. Les donnees m2 sont probablement inexploitables en l'etat (trop peu de reads / trop de bruit)." >&2
fi

########################
# 4. m3 - deja present, on ne relance pas DADA2
########################
echo "=== m3 : DADA2 deja present, aucune action ==="
if [[ ! -f "$WORKDIR/dada2/table_m3.qza" || ! -f "$WORKDIR/dada2/rep-seqs_m3.qza" ]]; then
  echo "ERREUR: m3 devrait deja avoir table/rep-seqs, verifier l'etat du dossier." >&2
fi

########################
# 5. Taxonomie m1 : diagnostic puis seuils relaches
########################
echo "=== m1 : diagnostic taxonomie (pourquoi tout Unassigned) ==="

REP_M1="$WORKDIR/dada2/rep-seqs_m1.qza"
REF_SEQS_M1="$WORKDIR/taxonomy/ref-seqs_m1.qza"
REF_TAX_M1="$WORKDIR/taxonomy/ref-taxonomy_m1.qza"
SEARCH_M1_OLD="$WORKDIR/taxonomy/search_m1.qza"

if [[ -f "$SEARCH_M1_OLD" ]]; then
  mkdir -p "$WORKDIR/tmp/search_m1_diag"
  qiime tools export \
    --input-path "$SEARCH_M1_OLD" \
    --output-path "$WORKDIR/tmp/search_m1_diag" \
    || true
  if [[ -f "$WORKDIR/tmp/search_m1_diag/blast6.tsv" ]]; then
    NB_HITS=$(wc -l < "$WORKDIR/tmp/search_m1_diag/blast6.tsv")
    echo "Nombre de lignes dans search_m1 (hits BLAST6 bruts) : $NB_HITS"
    echo "Si ce nombre est > 0, les sequences ont bien des hits, mais le seuil de consensus/identite les a rejetes."
    echo "Si ce nombre est proche de 0, la base de reference pla est probablement trop restreinte ou mal formatee."
  fi
fi

echo "=== m1 : relance taxonomie avec seuils relaches ==="

[[ -f "$WORKDIR/taxonomy/taxonomy_m1.qza" ]] && mv "$WORKDIR/taxonomy/taxonomy_m1.qza" "$WORKDIR/taxonomy/taxonomy_m1.qza.bak_$(date +%Y%m%d_%H%M%S)"
[[ -f "$WORKDIR/taxonomy/search_m1.qza" ]] && mv "$WORKDIR/taxonomy/search_m1.qza" "$WORKDIR/taxonomy/search_m1.qza.bak_$(date +%Y%m%d_%H%M%S)"

if [[ -f "$REP_M1" && -f "$REF_SEQS_M1" && -f "$REF_TAX_M1" ]]; then
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "$REP_M1" \
    --i-reference-reads "$REF_SEQS_M1" \
    --i-reference-taxonomy "$REF_TAX_M1" \
    --p-perc-identity 0.80 \
    --p-min-consensus 0.51 \
    --p-top-hits-only \
    --p-maxaccepts 10 \
    --p-threads "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m1.qza" \
    --o-search-results "$WORKDIR/taxonomy/search_m1.qza" \
    --verbose \
    2>&1 | tee "$WORKDIR/logs/taxonomy_m1_relaxed.log"
else
  echo "ATTENTION: taxonomie m1 non relancee (fichiers manquants)." >&2
fi

########################
# 6. Taxonomie m2 (si DADA2 m2 a produit des rep-seqs)
########################
echo "=== m2 : taxonomie si rep-seqs disponibles ==="

REP_M2="$WORKDIR/dada2/rep-seqs_m2.qza"
REF_SEQS_M2="$WORKDIR/taxonomy/ref-seqs_m2.qza"
REF_TAX_M2="$WORKDIR/taxonomy/ref-taxonomy_m2.qza"

if [[ -f "$REP_M2" && -f "$REF_SEQS_M2" && -f "$REF_TAX_M2" ]]; then
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "$REP_M2" \
    --i-reference-reads "$REF_SEQS_M2" \
    --i-reference-taxonomy "$REF_TAX_M2" \
    --p-perc-identity 0.80 \
    --p-min-consensus 0.51 \
    --p-top-hits-only \
    --p-maxaccepts 10 \
    --p-threads "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m2.qza" \
    --o-search-results "$WORKDIR/taxonomy/search_m2.qza" \
    --verbose \
    2>&1 | tee "$WORKDIR/logs/taxonomy_m2_relaxed.log"
else
  echo "ATTENTION: taxonomie m2 non relancee (rep-seqs_m2 absent ou bases manquantes)." >&2
fi

########################
# 7. Taxonomie m3 : chemin du classificateur SILVA CORRIGE
########################
echo "=== m3 : taxonomie avec chemin classificateur corrige ==="

# Le classificateur est reellement ici, pas dans WORKDIR/taxonomy :
M3_CLASSIFIER="$TAXODIR_GLOBAL/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
REP_M3="$WORKDIR/dada2/rep-seqs_m3.qza"

if [[ -f "$M3_CLASSIFIER" && -f "$REP_M3" ]]; then
  qiime feature-classifier classify-sklearn \
    --i-classifier "$M3_CLASSIFIER" \
    --i-reads "$REP_M3" \
    --p-n-jobs "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m3.qza" \
    --verbose \
    2>&1 | tee "$WORKDIR/logs/taxonomy_m3_classify.log"
else
  echo "ATTENTION : taxonomie m3 non faite. Verifier :" >&2
  echo "  - classificateur attendu : $M3_CLASSIFIER (existe: $([[ -f "$M3_CLASSIFIER" ]] && echo oui || echo NON))" >&2
  echo "  - rep-seqs_m3 attendu    : $REP_M3 (existe: $([[ -f "$REP_M3" ]] && echo oui || echo NON))" >&2
fi

########################
# 8. Controles negatifs : inspection visuelle
########################
echo "=== Inspection des controles ==="

for marker in m1 m2 m3; do
  TABLE="$WORKDIR/dada2/table_${marker}.qza"
  TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
  CTRL_TABLE="$WORKDIR/decontam/controls_table_${marker}.qza"
  CTRL_BAR="$WORKDIR/decontam/controls_barplot_${marker}.qzv"

  [[ -f "$TABLE" ]] || { echo "m${marker#m} : pas de table, controle ignore." ; continue; }
  [[ -f "$TAXO" ]] || { echo "$marker : pas de taxonomie, barplot controle ignore." ; continue; }

  qiime feature-table filter-samples \
    --i-table "$TABLE" \
    --m-metadata-file "$META" \
    --p-where "[sample-or-control]='control' AND [marker]='${marker}'" \
    --o-filtered-table "$CTRL_TABLE" \
    2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true

  if [[ -f "$CTRL_TABLE" ]]; then
    qiime taxa barplot \
      --i-table "$CTRL_TABLE" \
      --i-taxonomy "$TAXO" \
      --m-metadata-file "$META" \
      --o-visualization "$CTRL_BAR" \
      2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true
  else
    echo "$marker : table de controles vide ou filtree a zero (verifier colonnes 'sample-or-control' et 'marker' dans $META)."
  fi
done

########################
# 9. Decontam : identification
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
    2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true

  if [[ -f "$SCORE" ]]; then
    if [[ -f "$REP" ]]; then
      qiime quality-control decontam-score-viz \
        --i-decontam-scores "$SCORE" \
        --i-table "$TABLE" \
        --i-rep-seqs "$REP" \
        --p-threshold 0.1 \
        --o-visualization "$SCORE_VIZ" \
        2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true
    else
      qiime quality-control decontam-score-viz \
        --i-decontam-scores "$SCORE" \
        --i-table "$TABLE" \
        --p-threshold 0.1 \
        --o-visualization "$SCORE_VIZ" \
        2>&1 | tee -a "$WORKDIR/logs/decontam_${marker}.log" || true
    fi
  else
    echo "$marker : decontam-identify n'a pas produit de score (pas assez de controles ou table trop petite)."
  fi
done

########################
# 10. Retrait des contaminants (sans decontam-remove, absent du plugin)
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
  [[ -f "$SCORE" ]] || { echo "$marker : pas de score decontam, table finale = table brute filtree sur echantillons uniquement." ; }

  if [[ -f "$SCORE" ]]; then
    qiime feature-table filter-features \
      --i-table "$TABLE" \
      --m-metadata-file "$SCORE" \
      --p-where '[p] > 0.1 OR [p] IS NULL' \
      --o-filtered-table "$TABLE_DC" \
      2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true

    if [[ -f "$TABLE_DC" ]]; then
      qiime feature-table filter-seqs \
        --i-data "$REP" \
        --i-table "$TABLE_DC" \
        --o-filtered-data "$REP_DC" \
        2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
    fi
  else
    cp "$TABLE" "$TABLE_DC"
    cp "$REP" "$REP_DC"
  fi

  if [[ -f "$TABLE_DC" ]]; then
    qiime feature-table filter-samples \
      --i-table "$TABLE_DC" \
      --m-metadata-file "$META" \
      --p-where "[sample-or-control]='sample' AND [marker]='${marker}'" \
      --o-filtered-table "$TABLE_FINAL" \
      2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
  fi

  if [[ -f "$TABLE_FINAL" && -f "$REP_DC" ]]; then
    qiime feature-table filter-seqs \
      --i-data "$REP_DC" \
      --i-table "$TABLE_FINAL" \
      --o-filtered-data "$REP_FINAL" \
      2>&1 | tee -a "$WORKDIR/logs/filter_${marker}.log" || true
  fi
done

########################
# 11. Exports finaux
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

  TABLE_TO_USE="$TABLE_RAW"
  REP_TO_USE="$REP_RAW"

  if [[ -f "$TABLE_FINAL" && -f "$REP_FINAL" ]]; then
    TABLE_TO_USE="$TABLE_FINAL"
    REP_TO_USE="$REP_FINAL"
  fi

  [[ -f "$TABLE_TO_USE" ]] || { echo "$marker : aucune table disponible, export ignore." ; continue; }
  [[ -f "$REP_TO_USE" ]] || { echo "$marker : aucune rep-seqs disponible, export ignore." ; continue; }

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

echo "=== Pipeline v3 termine ==="
