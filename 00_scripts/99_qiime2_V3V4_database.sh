#!/bin/bash
#SBATCH --job-name=V3V4
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --cpus-per-task=36
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=/home/plstenge/Grand_Saint_Antoine/00_scripts/99_qiime2_V3V4_database.err
#SBATCH --output=/home/plstenge/Grand_Saint_Antoine/00_scripts/99_qiime2_V3V4_database.out

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
module load conda/4.12.0
source ~/.bashrc
conda activate rachis-qiime2-2026.7

cd /home/plstenge/Grand_Saint_Antoine/taxonomy

# Code from: https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494
qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

# If you'd like to be able to jump to steps that only take FeatureData[Sequence] as input you can convert your data to FeatureData[Sequence] like so:
qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

# “Culling” low-quality sequences with cull-seqs
# Here we’ll remove sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) and any homopolymers that are 8 or more bases in length. These are the default parameters.
qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

# Filtering sequences by length and taxonomy
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza 

#Dereplicating in uniq mode
qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.2-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza


qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
  --i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.2-ssu-nr99-classifier.qza    
