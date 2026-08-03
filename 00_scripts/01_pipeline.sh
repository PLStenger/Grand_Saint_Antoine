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
# Pipeline QIIME2 - Analyse multi-marqueurs (pla, caf, V4V5) avec gestion
#                    des contrôles négatifs (NTC / Neg)
# Environnement : conda activate rachis-qiime2-2026.7
#############################################################################

########################
# 0. Paramètres généraux
########################
BASEDIR="/storage/groups/gdec/shared_paleo/filtermask"
WORKDIR="/home/plstenge/Grand_Saint_Antoine"
THREADS=16

mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata}
cd "$WORKDIR"

# Marqueurs et primers (truseq stub retiré : on garde uniquement la partie
# spécifique du primer après les Ns, cutadapt gère le reste via --p-front)
declare -A PRIMER_F
declare -A PRIMER_R

PRIMER_F[m1]="GACTGGGTTCGGGCACATG"          # pla   F
PRIMER_R[m1]="AGACTTTGGCATTAGGTGTG"         # pla   R
PRIMER_F[m2]="aaccagcccgcatcactctta"        # caf1  F
PRIMER_R[m2]="atcacccgcggcatctgta"          # caf1  R
PRIMER_F[m3]="GTGYCAGCMGCCGCGGTAA"          # V4V5  F (515F-Y)
PRIMER_R[m3]="CCGYCAATTYMTTTRAGTTT"         # V4V5  R (926R)

MARKERS=(m1 m2 m3)

#####################################################
# 1. Construction des manifests (un par marqueur)
#####################################################
# Structure attendue : $BASEDIR/<ID>-<Rep>-<marker>/<sampleid>_R1.fastq.gz
# ex: $BASEDIR/1A-m1/1a-m1_R1.fastq.gz  (les noms internes varient un peu,
# ex ntc6-m7 -> on cherche donc dynamiquement les fastq dans chaque dossier)

for marker in "${MARKERS[@]}"; do
  manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
  echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$manifest"

  for d in "$BASEDIR"/*-"$marker"; do
    [ -d "$d" ] || continue
    sampledir=$(basename "$d")                 # ex: 1A-m1, NTC-A-m1, Neg-m1
    r1=$(find "$d" -maxdepth 1 -iname "*_R1*.fastq.gz" | head -n1)
    r2=$(find "$d" -maxdepth 1 -iname "*_R2*.fastq.gz" | head -n1)
    if [[ -z "$r1" || -z "$r2" ]]; then
      echo "ATTENTION: fastq manquants dans $d" >&2
      continue
    fi
    # sample-id QIIME2 = nom du dossier (garanti unique)
    echo -e "${sampledir}\t${r1}\t${r2}" >> "$manifest"
  done
  echo "Manifest ${marker} : $manifest ($(($(wc -l < "$manifest")-1)) échantillons)"
done

#####################################################
# 2. Fichier de métadonnées (commun, avec statut contrôle)
#####################################################
meta="$WORKDIR/metadata/sample-metadata.tsv"
echo -e "sample-id\tsample\treplicate\tmarker\tsample-or-control" > "$meta"

for marker in "${MARKERS[@]}"; do
  tail -n +2 "$WORKDIR/manifests/manifest_${marker}.tsv" | cut -f1 | while read -r sid; do
    if [[ "$sid" =~ ^(NTC|Neg) ]]; then
      status="control"
    else
      status="sample"
    fi
    # extraction sample (chiffre) et replicat (lettre) quand applicable
    samplenum=$(echo "$sid" | grep -oP '^[0-9]+' || echo "NA")
    repl=$(echo "$sid" | grep -oP '(?<=[0-9])[A-C]|(?<=NTC-)[A-C]' || echo "NA")
    echo -e "${sid}\t${samplenum}\t${repl}\t${marker}\t${status}" >> "$meta"
  done
done
echo "Metadata écrite : $meta"

#####################################################
# 3. Import QIIME2 + suppression primers (cutadapt) + DADA2, par marqueur
#####################################################
for marker in "${MARKERS[@]}"; do
  echo "=== Traitement marqueur $marker ==="
  manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
  demux="$WORKDIR/imported/demux_${marker}.qza"
  trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
  trimmed_qzv="$WORKDIR/trimmed/trimmed_${marker}.qzv"

  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$manifest" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$demux"

  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "$demux" \
    --p-cores "$THREADS" \
    --p-front-f "${PRIMER_F[$marker]}" \
    --p-front-r "${PRIMER_R[$marker]}" \
    --p-match-read-wildcards \
    --p-discard-untrimmed \
    --p-no-indels \
    --o-trimmed-sequences "$trimmed" \
    --verbose > "$WORKDIR/trimmed/cutadapt_${marker}.log"

  qiime demux summarize \
    --i-data "$trimmed" \
    --o-visualization "$trimmed_qzv"

  echo "-> Inspecter $trimmed_qzv pour choisir les longueurs de troncature DADA2"
done

#####################################################
# 4. DADA2 (denoise-paired) par marqueur
#    IMPORTANT : ouvrez les .qzv de l'étape 3 (qiime tools view) pour
#    ajuster --p-trunc-len-f/--p-trunc-len-r selon la qualité + la longueur
#    attendue de l'amplicon (pla/caf ~ courts, V4V5 ~ 400 pb).
#    Les valeurs ci-dessous sont des points de départ à valider.
#####################################################
declare -A TRUNC_F
declare -A TRUNC_R
TRUNC_F[m1]=0   # pla   - ajuster après inspection qzv (0 = pas de troncature)
TRUNC_R[m1]=0
TRUNC_F[m2]=0   # caf
TRUNC_R[m2]=0
TRUNC_F[m3]=280 # V4V5
TRUNC_R[m3]=220

for marker in "${MARKERS[@]}"; do
  trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
  table="$WORKDIR/dada2/table_${marker}.qza"
  repseqs="$WORKDIR/dada2/rep-seqs_${marker}.qza"
  stats="$WORKDIR/dada2/stats_${marker}.qza"

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
#    - V4V5 (m3) : classificateur SILVA (à télécharger au préalable, ex.
#      silva-138-99-nb-classifier.qza dans $WORKDIR/taxonomy/)
#    - pla / caf (m1, m2) : marqueurs spécifiques Yersinia, pas de base
#      SILVA adaptée -> utiliser classify-consensus-blast/vsearch contre une
#      base de référence dédiée (ex. séquences pla/caf1 de Y. pestis/
#      Enterobacteriaceae) que vous devez fournir.
#####################################################

# --- V4V5 : classification SILVA ---
CLASSIFIER="$WORKDIR/taxonomy/silva-138-99-nb-classifier.qza"
if [[ -f "$CLASSIFIER" ]]; then
  qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$WORKDIR/dada2/rep-seqs_m3.qza" \
    --p-n-jobs "$THREADS" \
    --o-classification "$WORKDIR/taxonomy/taxonomy_m3.qza"
else
  echo "ATTENTION : classificateur SILVA absent, téléchargez-le avant l'étape taxonomie V4V5"
fi

# --- pla / caf : classification via base de référence dédiée (à adapter) ---
for marker in m1 m2; do
  REF_SEQS="$WORKDIR/taxonomy/ref-seqs_${marker}.qza"      # FeatureData[Sequence]
  REF_TAX="$WORKDIR/taxonomy/ref-taxonomy_${marker}.qza"   # FeatureData[Taxonomy]
  if [[ -f "$REF_SEQS" && -f "$REF_TAX" ]]; then
    qiime feature-classifier classify-consensus-vsearch \
      --i-query "$WORKDIR/dada2/rep-seqs_${marker}.qza" \
      --i-reference-reads "$REF_SEQS" \
      --i-reference-taxonomy "$REF_TAX" \
      --p-threads "$THREADS" \
      --o-classification "$WORKDIR/taxonomy/taxonomy_${marker}.qza" \
      --o-search-results "$WORKDIR/taxonomy/search_${marker}.qza"
  else
    echo "ATTENTION : fournissez ref-seqs_${marker}.qza + ref-taxonomy_${marker}.qza (séquences pla/caf1 annotées) pour classer $marker"
  fi
done

#####################################################
# 6. Inspection des contrôles AVANT nettoyage
#    -> barplot taxonomique restreint aux échantillons contrôles
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

  echo "-> Ouvrez controls_barplot_${marker}.qzv (qiime tools view) pour voir ce qu'il y a dans les NTC/Neg de $marker"
done

#####################################################
# 7. Décontamination (decontam, méthode "prevalence" via metadata)
#    -> identifie les ASV significativement associées aux contrôles
#    -> les retire des échantillons réels
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

  # Retrait des features identifiées comme contaminantes (score <= 0.1)
  qiime quality-control decontam-remove \
    --i-decontam-scores "$WORKDIR/decontam/decontam-scores_${marker}.qza" \
    --i-table "$table" \
    --i-rep-seqs "$repseqs" \
    --p-threshold 0.1 \
    --o-filtered-table "$WORKDIR/decontam/table_${marker}_decontam.qza" \
    --o-filtered-rep-seqs "$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"

  echo "-> Inspectez decontam-scoreviz_${marker}.qzv avant de valider le seuil 0.1 (ajustable)"
done

#####################################################
# 8. Retrait des échantillons contrôles eux-mêmes de la table finale
#    (une fois les contaminants ASV retirés des vrais échantillons)
#####################################################
for marker in "${MARKERS[@]}"; do
  table_dc="$WORKDIR/decontam/table_${marker}_decontam.qza"
  [ -f "$table_dc" ] || continue

  qiime feature-table filter-samples \
    --i-table "$table_dc" \
    --m-metadata-file "$meta" \
    --p-where "[marker]='${marker}' AND [sample-or-control]='sample'" \
    --o-filtered-table "$WORKDIR/decontam/table_${marker}_final.qza"

  qiime feature-table filter-seqs \
    --i-data "$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza" \
    --i-table "$WORKDIR/decontam/table_${marker}_final.qza" \
    --o-filtered-data "$WORKDIR/decontam/rep-seqs_${marker}_final.qza"
done

#####################################################
# 9. Table taxonomique finale + exports (pas de raréfaction, comme demandé)
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
  qiime tools export --input-path "$taxo"       --output-path "$WORKDIR/exports/${marker}_final"

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

  echo "=== Marqueur ${marker} : table ASV finale + taxonomie exportées dans exports/${marker}_final/ ==="
done

echo "Pipeline terminé."
