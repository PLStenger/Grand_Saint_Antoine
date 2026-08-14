#!/bin/bash
#SBATCH --job-name=GSA
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline_full.err
#SBATCH --output=/home/plstenge/Grand_Saint_Antoine/00_scripts/01_pipeline_full.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

#############################################################################
# FULL PIPELINE - Grand Saint Antoine - multi-marker analysis (m1=pla, m2=caf1,
# m3=V4V5), with negative control handling (NTC / Neg).
#
# This single script merges and fixes the two previous scripts:
#   - 01_pipeline.sh          (initial pipeline, worked up to a point)
#   - 01_pipeline_resume-2.sh (patched re-run: fixed m2 DADA2 strategy,
#                              fixed SILVA classifier path for m3, relaxed
#                              taxonomy thresholds for m1/m2, manual
#                              contaminant filtering since "decontam-remove"
#                              does not exist as a q2 action)
#
# NEW in this version (per user request):
#   - Illumina universal adapter file auto-generated (ADAPTERFILE)
#   - Trimmomatic quality/adapter trimming step BEFORE QIIME2 import
#   - FastQC on raw reads + MultiQC aggregated PER MARKER (m1, m2, m3)
#   - FastQC on trimmomatic-cleaned reads + MultiQC aggregated PER MARKER
#############################################################################

########################
# 0. General parameters
########################
#BASEDIR="/storage/groups/gdec/shared_paleo/E1739/filtermask"
BASEDIR="/storage/groups/gdec/shared_paleo/E1739-Ps12-2/DefaultProject"
#WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_filtermask_pipeline"
WORKDIR="/home/plstenge/Grand_Saint_Antoine/02_E1739_Ps12_2_DefaultProject"
# Global taxonomy directory: some reference files (e.g. the SILVA classifier)
# actually live here rather than inside WORKDIR/taxonomy (bug found during
# the first run and fixed in the resume script).
TAXODIR_GLOBAL="/home/plstenge/Grand_Saint_Antoine/taxonomy"
THREADS=16
TRIM_THREADS=8
JAVA_MEM="60G"

# Set to 1 to also strip generic Illumina adapters during the QIIME2 cutadapt
# primer-removal step (kept OFF by default, exactly like the original
# pipeline, because it previously caused an empty-adapter pass-through bug
# in q2-cutadapt). Trimmomatic (below) already handles adapter removal on
# raw reads, so this is mostly redundant safety.
USE_ILLUMINA_ADAPTER_TRIM=0

# Create the full working directory tree, including the new raw/cleaned
# QC folders and the trimmed-fastq folder used by Trimmomatic.
mkdir -p "$WORKDIR"/{manifests,imported,trimmed,dada2,taxonomy,decontam,exports,metadata,logs,tmp}
mkdir -p "$WORKDIR"/qc/{raw,cleaned}/{m1,m2,m3}
mkdir -p "$WORKDIR"/trimmomatic
cd "$WORKDIR"

META="$WORKDIR/metadata/sample-metadata.tsv"

# Primers specific to each marker (used later by cutadapt for primer removal)
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

# Generic Illumina adapters (kept for optional use in q2-cutadapt step)
ILLUMINA_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ILLUMINA_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

#####################################################
# 0bis. Sample directory list (explicit, GSA project)
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

#####################################################
# 1. Build the Illumina universal adapter FASTA file
#    (used by Trimmomatic ILLUMINACLIP)
#####################################################
ADAPTERFILE="$WORKDIR/trimmomatic/illumina_universal_adapters.fa"

echo "=== Writing Illumina universal adapter file: $ADAPTERFILE ==="
cat > "$ADAPTERFILE" << 'ADAPTER_EOF'
>PrefixPE/1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>PrefixPE/2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_fwd
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_rev
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Nextera_Trans1
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Nextera_Trans2
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>NexteraPE-PE/1
CTGTCTCTTATACACATCT
>NexteraPE-PE/2
CTGTCTCTTATACACATCT
>Illumina_Single_End_Adapter_1
GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>Illumina_Single_End_Adapter_2
GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>Illumina_Paired_End_Adapter_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Paired_End_Adapter_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_Multiplexing_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Multiplexing_Adapter_2
GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_Multiplexing_Index_Sequencing_Primer
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Illumina_Multiplexing_Read2_Sequencing_Primer
GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_DpnII_expression_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_DpnII_expression_Adapter_2
CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
>Illumina_DpnII_Gex_Sequencing_Primer
CGACAGGTTCAGAGTTCTACAGTCCGACGATC
>Illumina_NlaIII_expression_Adapter_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_NlaIII_expression_Adapter_2
CATGCAGAAGACGGCATACGAGCTCTTCCGATCT
>Illumina_NlaIII_Gex_Sequencing_Primer
CCGACTATGCCGTCTGTTCCGAAGGTCCGACGATC
>Illumina_Small_RNA_Adapter_1
ATCTCGTATGCCGTCTTCTGCTTG
>Illumina_Small_RNA_RT_Primer
CAAGCAGAAGACGGCATACGA
>Illumina_Small_RNA_PCR_Primer_1
CAAGCAGAAGACGGCATACGA
>Illumina_Small_RNA_PCR_Primer_2
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ADAPTER_EOF

[[ -s "$ADAPTERFILE" ]] || { echo "ERROR: adapter file was not created correctly." >&2; exit 1; }
echo "Adapter file ready with $(grep -c '^>' "$ADAPTERFILE") sequences."

#####################################################
# 2. Resolve raw R1/R2 fastq paths for every sample
#    (shared helper used by QC steps, Trimmomatic and
#    later by the manifest-building step)
#####################################################
declare -A RAW_R1
declare -A RAW_R2

echo "=== Resolving raw fastq paths for all samples ==="
for sampledir in "${GSA_DIRS[@]}"; do
    d="$BASEDIR/$sampledir"
    if [[ ! -d "$d" ]]; then
        echo "WARNING: directory missing: $d" >&2
        continue
    fi
    r1=$(find "$d" -maxdepth 1 -type f -iname "*_R1*.fastq.gz" | sort | head -n1)
    r2=$(find "$d" -maxdepth 1 -type f -iname "*_R2*.fastq.gz" | sort | head -n1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: fastq files missing in $d" >&2
        continue
    fi
    RAW_R1[$sampledir]="$r1"
    RAW_R2[$sampledir]="$r2"
done

#####################################################
# 3. FastQC + MultiQC on RAW data, PER MARKER
#####################################################
echo "=== FastQC on raw reads ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/raw/${marker}"
    mkdir -p "$outdir"
    for sampledir in "${GSA_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${RAW_R1[$sampledir]:-}" && -n "${RAW_R2[$sampledir]:-}" ]] || continue
        fastqc -t "$TRIM_THREADS" -o "$outdir" "${RAW_R1[$sampledir]}" "${RAW_R2[$sampledir]}" \
            2>&1 | tee -a "$WORKDIR/logs/fastqc_raw_${marker}.log"
    done
done

echo "=== MultiQC on raw reads, one report per marker ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/raw/${marker}"
    multiqc "$outdir" -o "$outdir" -n "multiqc_raw_${marker}" \
        2>&1 | tee -a "$WORKDIR/logs/multiqc_raw_${marker}.log"
done

#####################################################
# 4. Trimmomatic cleaning (quality trimming + adapter
#    removal) for every sample, using the generated
#    ADAPTERFILE
#####################################################
echo "=== Trimmomatic cleaning ==="
declare -A CLEAN_R1
declare -A CLEAN_R2

for sampledir in "${GSA_DIRS[@]}"; do
    [[ -n "${RAW_R1[$sampledir]:-}" && -n "${RAW_R2[$sampledir]:-}" ]] || continue

    out_paired_r1="$WORKDIR/trimmomatic/${sampledir}_R1.paired.fastq.gz"
    out_single_r1="$WORKDIR/trimmomatic/${sampledir}_R1.single.fastq.gz"
    out_paired_r2="$WORKDIR/trimmomatic/${sampledir}_R2.paired.fastq.gz"
    out_single_r2="$WORKDIR/trimmomatic/${sampledir}_R2.single.fastq.gz"

    # Parameters adapted from the user's example script:
    # ILLUMINACLIP:ADAPTERFILE:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:26:30 MINLEN:150
    trimmomatic PE -Xmx"$JAVA_MEM" -threads "$TRIM_THREADS" -phred33 \
        "${RAW_R1[$sampledir]}" "${RAW_R2[$sampledir]}" \
        "$out_paired_r1" "$out_single_r1" \
        "$out_paired_r2" "$out_single_r2" \
        ILLUMINACLIP:"$ADAPTERFILE":2:30:10 \
        LEADING:30 TRAILING:30 SLIDINGWINDOW:26:30 MINLEN:150 \
        2>&1 | tee "$WORKDIR/logs/trimmomatic_${sampledir}.log"

    if [[ -s "$out_paired_r1" && -s "$out_paired_r2" ]]; then
        CLEAN_R1[$sampledir]="$out_paired_r1"
        CLEAN_R2[$sampledir]="$out_paired_r2"
    else
        echo "ERROR: Trimmomatic did not produce paired output for $sampledir" >&2
    fi
done

#####################################################
# 5. FastQC + MultiQC on CLEANED data, PER MARKER
#####################################################
echo "=== FastQC on Trimmomatic-cleaned reads ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/cleaned/${marker}"
    mkdir -p "$outdir"
    for sampledir in "${GSA_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${CLEAN_R1[$sampledir]:-}" && -n "${CLEAN_R2[$sampledir]:-}" ]] || continue
        fastqc -t "$TRIM_THREADS" -o "$outdir" "${CLEAN_R1[$sampledir]}" "${CLEAN_R2[$sampledir]}" \
            2>&1 | tee -a "$WORKDIR/logs/fastqc_cleaned_${marker}.log"
    done
done

echo "=== MultiQC on cleaned reads, one report per marker ==="
for marker in "${MARKERS[@]}"; do
    outdir="$WORKDIR/qc/cleaned/${marker}"
    multiqc "$outdir" -o "$outdir" -n "multiqc_cleaned_${marker}" \
        2>&1 | tee -a "$WORKDIR/logs/multiqc_cleaned_${marker}.log"
done

#####################################################
# 6. Build QIIME2 manifests from the CLEANED
#    (Trimmomatic-paired) fastq files
#####################################################
echo "=== Building manifests from cleaned reads ==="
for marker in "${MARKERS[@]}"; do
    manifest="$WORKDIR/manifests/manifest_${marker}.tsv"
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$manifest"

    for sampledir in "${GSA_DIRS[@]}"; do
        [[ "$sampledir" == *-"$marker" ]] || continue
        [[ -n "${CLEAN_R1[$sampledir]:-}" && -n "${CLEAN_R2[$sampledir]:-}" ]] || continue
        printf "%s\t%s\t%s\n" "$sampledir" "${CLEAN_R1[$sampledir]}" "${CLEAN_R2[$sampledir]}" >> "$manifest"
    done

    echo "Manifest ${marker}: $manifest"
    column -t -s $'\t' "$manifest" | head
    echo
done

#####################################################
# 7. Metadata file (shared across markers), including
#    control status
#####################################################
printf "sample-id\tsample\treplicate\tmarker\tsample-or-control\n" > "$META"

for sampledir in "${GSA_DIRS[@]}"; do
    # Split GSA directory name, e.g. "1A-m3"
    sample="${sampledir%%-*}"      # 1A
    marker="${sampledir##*-}"      # m3
    replicate="${sample: -1}"      # A
    core="${sample::-1}"           # 1
    soc="sample"

    case "$sampledir" in
        NTC-*|Neg-*)
            soc="control"
            ;;
    esac

    printf "%s\t%s\t%s\t%s\t%s\n" "$sampledir" "$core" "$replicate" "$marker" "$soc" >> "$META"
done

echo "Metadata written: $META"

#####################################################
# 8. QIIME2 import + primer/adapter removal (cutadapt)
#    NOTE: reads are already Trimmomatic-cleaned, so
#    this step only removes the marker-specific primers.
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

    [[ -f "$trimmed" ]] || { echo "ERROR: cutadapt did not create $trimmed" >&2; exit 1; }

    qiime demux summarize \
        --i-data "$trimmed" \
        --o-visualization "$trimmed_qzv"
done

#####################################################
# 9. DADA2 denoise-paired, per marker
#    m2 uses the alternative strategy validated in the
#    resume script (pseudo pooling, retry with/without
#    chimera detection) because the default independent
#    pooling produced "No features remain after denoising".
#####################################################
declare -A TRUNC_F
declare -A TRUNC_R
declare -A MAXEE_F
declare -A MAXEE_R

TRUNC_F[m1]=0
TRUNC_R[m1]=0
MAXEE_F[m1]=5.0
MAXEE_R[m1]=5.0

TRUNC_F[m2]=0
TRUNC_R[m2]=0
MAXEE_F[m2]=10.0
MAXEE_R[m2]=10.0

TRUNC_F[m3]=220
TRUNC_R[m3]=180
MAXEE_F[m3]=5.0
MAXEE_R[m3]=5.0

# --- m1 and m3: standard DADA2 run ---
for marker in m1 m3; do
    trimmed="$WORKDIR/trimmed/trimmed_${marker}.qza"
    table="$WORKDIR/dada2/table_${marker}.qza"
    repseqs="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    stats="$WORKDIR/dada2/stats_${marker}.qza"
    bt_stats="$WORKDIR/dada2/base-transition-stats_${marker}.qza"

    [[ -f "$trimmed" ]] || { echo "ERROR: missing file $trimmed" >&2; continue; }

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$trimmed" \
        --p-trunc-len-f "${TRUNC_F[$marker]}" \
        --p-trunc-len-r "${TRUNC_R[$marker]}" \
        --p-max-ee-f "${MAXEE_F[$marker]}" \
        --p-max-ee-r "${MAXEE_R[$marker]}" \
        --p-n-threads "$THREADS" \
        --o-table "$table" \
        --o-representative-sequences "$repseqs" \
        --o-denoising-stats "$stats" \
        --o-base-transition-stats "$bt_stats" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/dada2_${marker}.log"

    if [[ -f "$table" ]]; then
        qiime feature-table summarize \
            --i-table "$table" \
            --m-metadata-file "$META" \
            --o-feature-frequencies "$WORKDIR/dada2/feature-freq_${marker}.qza" \
            --o-sample-frequencies "$WORKDIR/dada2/sample-freq_${marker}.qza" \
            --o-summary "$WORKDIR/dada2/table_${marker}.qzv"

        qiime feature-table tabulate-seqs \
            --i-data "$repseqs" \
            --o-visualization "$WORKDIR/dada2/rep-seqs_${marker}.qzv"
    else
        echo "ERROR: DADA2 did not produce a table for $marker" >&2
    fi
done

# --- m2: alternative strategy (pseudo pooling) ---
echo "=== m2: DADA2 with pseudo pooling (attempt 1: no chimera removal) ==="
M2_TRUNC_F=0
M2_TRUNC_R=0
M2_MAXEE_F=10.0
M2_MAXEE_R=10.0

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$WORKDIR/trimmed/trimmed_m2.qza" \
    --p-trunc-len-f "$M2_TRUNC_F" \
    --p-trunc-len-r "$M2_TRUNC_R" \
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
    echo "Attempt 1 (pseudo, no chimera) failed, attempt 2: pseudo pooling + consensus chimera removal" >&2

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$WORKDIR/trimmed/trimmed_m2.qza" \
        --p-trunc-len-f "$M2_TRUNC_F" \
        --p-trunc-len-r "$M2_TRUNC_R" \
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
    echo "m2: DADA2 produced a table (check logs to see whether chimera removal was applied)."
else
    echo "ERROR: m2 DADA2 failed even with pseudo pooling. Data may be unusable (too few / too noisy reads)." >&2
fi

#####################################################
# 10. Taxonomy
#     - m3 (V4V5): SILVA sklearn classifier
#       (path checked both in WORKDIR/taxonomy and in
#       TAXODIR_GLOBAL, since the classifier was found
#       to actually live in the latter)
#     - m1 (pla) / m2 (caf1): custom NCBI-built reference
#       database + classify-consensus-vsearch, with
#       relaxed thresholds (perc-identity 0.80) since the
#       default 0.97 produced 100% Unassigned on m1
#####################################################

# --- m3: locate SILVA classifier, trying both known paths ---
CLASSIFIER_LOCAL="$WORKDIR/taxonomy/silva-138.2-ssu-nr99-515f-926r-classifier.qza"
CLASSIFIER_GLOBAL="$TAXODIR_GLOBAL/silva-138.2-ssu-nr99-515f-926r-classifier.qza"

if [[ -f "$CLASSIFIER_LOCAL" ]]; then
    M3_CLASSIFIER="$CLASSIFIER_LOCAL"
elif [[ -f "$CLASSIFIER_GLOBAL" ]]; then
    M3_CLASSIFIER="$CLASSIFIER_GLOBAL"
else
    M3_CLASSIFIER=""
fi

REP_M3="$WORKDIR/dada2/rep-seqs_m3.qza"

if [[ -n "$M3_CLASSIFIER" && -f "$REP_M3" ]]; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "$M3_CLASSIFIER" \
        --i-reads "$REP_M3" \
        --p-n-jobs "$THREADS" \
        --o-classification "$WORKDIR/taxonomy/taxonomy_m3.qza" \
        --verbose \
        2>&1 | tee "$WORKDIR/logs/taxonomy_m3.log"
else
    echo "WARNING: m3 taxonomy skipped. Classifier found: ${M3_CLASSIFIER:-NONE}; rep-seqs found: $([[ -f "$REP_M3" ]] && echo yes || echo no)" >&2
fi

# --- m1 / m2: build custom NCBI reference database if missing ---
declare -A NCBI_QUERY
NCBI_QUERY[m1]='(Yersinia pestis[Organism] AND pla[Gene]) AND (plasminogen activator[Title] OR pla[Title])'
NCBI_QUERY[m2]='(Yersinia pestis[Organism] OR Enterobacteriaceae[Organism]) AND (caf1[Gene] OR "capsular antigen F1"[Title] OR caf1[Title])'

for marker in m1 m2; do
    REF_SEQS="$WORKDIR/taxonomy/ref-seqs_${marker}.qza"
    REF_TAX="$WORKDIR/taxonomy/ref-taxonomy_${marker}.qza"

    if [[ ! -f "$REF_SEQS" || ! -f "$REF_TAX" ]]; then
        echo "Reference database for ${marker} missing -> building it from NCBI..."
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

    REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"

    if [[ -f "$REP" && -f "$REF_SEQS" && -f "$REF_TAX" ]]; then
        qiime feature-classifier classify-consensus-vsearch \
            --i-query "$REP" \
            --i-reference-reads "$REF_SEQS" \
            --i-reference-taxonomy "$REF_TAX" \
            --p-perc-identity 0.80 \
            --p-min-consensus 0.51 \
            --p-top-hits-only \
            --p-maxaccepts 10 \
            --p-threads "$THREADS" \
            --o-classification "$WORKDIR/taxonomy/taxonomy_${marker}.qza" \
            --o-search-results "$WORKDIR/taxonomy/search_${marker}.qza" \
            --verbose \
            2>&1 | tee "$WORKDIR/logs/taxonomy_${marker}.log"
    else
        echo "WARNING: ${marker} classification skipped (rep-seqs or reference database unavailable)." >&2
    fi
done

#####################################################
# 11. Inspect negative controls BEFORE decontamination
#####################################################
for marker in "${MARKERS[@]}"; do
    TABLE="$WORKDIR/dada2/table_${marker}.qza"
    TAXO="$WORKDIR/taxonomy/taxonomy_${marker}.qza"
    [[ -f "$TABLE" && -f "$TAXO" ]] || continue

    qiime feature-table filter-samples \
        --i-table "$TABLE" \
        --m-metadata-file "$META" \
        --p-where "[marker]='${marker}' AND [sample-or-control]='control'" \
        --o-filtered-table "$WORKDIR/decontam/controls_table_${marker}.qza" \
        2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true

    if [[ -f "$WORKDIR/decontam/controls_table_${marker}.qza" ]]; then
        qiime taxa barplot \
            --i-table "$WORKDIR/decontam/controls_table_${marker}.qza" \
            --i-taxonomy "$TAXO" \
            --m-metadata-file "$META" \
            --o-visualization "$WORKDIR/decontam/controls_barplot_${marker}.qzv" \
            2>&1 | tee -a "$WORKDIR/logs/controls_${marker}.log" || true
    fi
done

#####################################################
# 12. Decontam identification (prevalence method)
#####################################################
for marker in "${MARKERS[@]}"; do
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
        echo "$marker: decontam-identify did not produce scores (not enough controls or table too small)."
    fi
done

#####################################################
# 13. Remove contaminant features and control samples
#     NOTE: "qiime quality-control decontam-remove" does
#     not exist in this q2 plugin version, so contaminant
#     features are removed manually via
#     feature-table filter-features on the p-value score.
#####################################################
for marker in "${MARKERS[@]}"; do
    TABLE="$WORKDIR/dada2/table_${marker}.qza"
    REP="$WORKDIR/dada2/rep-seqs_${marker}.qza"
    SCORE="$WORKDIR/decontam/decontam-scores_${marker}.qza"

    TABLE_DC="$WORKDIR/decontam/table_${marker}_decontam.qza"
    REP_DC="$WORKDIR/decontam/rep-seqs_${marker}_decontam.qza"
    TABLE_FINAL="$WORKDIR/decontam/table_${marker}_final.qza"
    REP_FINAL="$WORKDIR/decontam/rep-seqs_${marker}_final.qza"

    [[ -f "$TABLE" && -f "$REP" ]] || continue

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
        echo "$marker: no decontam score available, using raw table/rep-seqs before sample filtering."
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

#####################################################
# 14. Final taxonomy barplots + exports
#     (falls back to the raw, non-decontaminated table
#     if the decontaminated final table is unavailable)
#####################################################
for marker in "${MARKERS[@]}"; do
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

    [[ -f "$TABLE_TO_USE" ]] || { echo "$marker: no table available, export skipped."; continue; }
    [[ -f "$REP_TO_USE" ]]   || { echo "$marker: no rep-seqs available, export skipped."; continue; }

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

    qiime tools export --input-path "$TABLE_TO_USE" --output-path "$EXPORT_DIR/table_export"
    qiime tools export --input-path "$REP_TO_USE" --output-path "$EXPORT_DIR/repseq_export"

    if [[ -f "$TAXO" ]]; then
        qiime tools export --input-path "$TAXO" --output-path "$EXPORT_DIR/taxonomy_export"
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

echo "=== Full pipeline finished ==="
