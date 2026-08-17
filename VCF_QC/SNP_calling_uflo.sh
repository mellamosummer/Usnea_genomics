#!/bin/bash
# ============================================================
#  SNP calling pipeline  (bwa mem -> markdup -> RG -> freebayes -> vcftools)
#  Expects a sample sheet:  /home/FM/sblanco/GWAS_8072026/samples.txt
#  Tab-separated, one line per sample:
#     <sample_id> <path/to/R1.fastq.gz> <path/to/R2.fastq.gz>
# ============================================================

WORKDIR="/home/FM/sblanco/GWAS_8072026"
REF_NAME="uflo"
MD_DIR="01_markdup"
RG_DIR="02_readgroups"
FB_DIR="03_freebayes"
FILT_DIR="04_filtering"

mkdir -p "${WORKDIR}/${MD_DIR}" \
         "${WORKDIR}/${RG_DIR}" \
         "${WORKDIR}/${FB_DIR}" \
         "${WORKDIR}/${FILT_DIR}/filter_steps" \
         "${WORKDIR}/logs"

# ---------- make BAM files with duplicates marked ----------

while IFS=$'\t' read -r sample R1 R2; do
  bwa-mem2 mem -t 20 -M "${WORKDIR}/reference/${REF_NAME}/genome.fa" \
    "${WORKDIR}/raw_reads/${R1}" "${WORKDIR}/raw_reads/${R2}" \
    | samtools sort -@ 20 -T "${WORKDIR}/${MD_DIR}/${sample}.sorttmp" \
      -o "${WORKDIR}/${MD_DIR}/${sample}.sorted.bam" -

  sambamba markdup --nthreads=20 --overflow-list-size 600000 \
    "${WORKDIR}/${MD_DIR}/${sample}.sorted.bam" \
    "${WORKDIR}/${MD_DIR}/${sample}.sorted.bam.md"

  samtools addreplacerg \
    -r "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
    -o "${WORKDIR}/${RG_DIR}/${sample}.rg.bam" \
    "${WORKDIR}/${MD_DIR}/${sample}.sorted.bam.md"

  samtools index "${WORKDIR}/${RG_DIR}/${sample}.rg.bam"
done < "${WORKDIR}/samples.txt"

# ---------- sanity check before variant calling ----------
samtools quickcheck -v "${WORKDIR}/${RG_DIR}"/*.rg.bam

# ---------- parallel freebayes ----------
fasta_generate_regions.py "${WORKDIR}/reference/${REF_NAME}/genome.fa.fai" 100000 \
  > "${WORKDIR}/${FB_DIR}/regions.txt"

freebayes-parallel "${WORKDIR}/${FB_DIR}/regions.txt" 16 \
  -f "${WORKDIR}/reference/${REF_NAME}/genome.fa" -p 1 \
  "${WORKDIR}/${RG_DIR}"/*.rg.bam \
  > "${WORKDIR}/${FB_DIR}/${REF_NAME}.vcf"

# ---------- VCFTools filtering ----------

# 6.1 biallelic SNPs only
vcftools --vcf "${WORKDIR}/${FB_DIR}/${REF_NAME}.vcf" \
  --remove-indels --min-alleles 2 --max-alleles 2 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/1_biallelic_snps.vcf.gz"

# 6.2 site quality
vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/1_biallelic_snps.vcf.gz" \
  --minQ 30 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/2_qual.vcf.gz"

# 6.3 depth (site mean + per-genotype)
vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/2_qual.vcf.gz" \
  --min-meanDP 10 --max-meanDP 250 \
  --minDP 10 --maxDP 250 \
  --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/3_depth.vcf.gz"

# 6.4 minor allele frequency
vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/3_depth.vcf.gz" \
  --maf 0.10 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz"

# 6.5 max-missing (last), one final VCF per threshold
vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz" \
  --max-missing 0.95 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_95.vcf.gz"
cp "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_95.vcf.gz" \
   "${WORKDIR}/${FILT_DIR}/${REF_NAME}_filtered_miss0_95.vcf.gz"

vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz" \
  --max-missing 0.90 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_90.vcf.gz"
cp "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_90.vcf.gz" \
   "${WORKDIR}/${FILT_DIR}/${REF_NAME}_filtered_miss0_90.vcf.gz"

vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz" \
  --max-missing 0.85 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_85.vcf.gz"
cp "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_85.vcf.gz" \
   "${WORKDIR}/${FILT_DIR}/${REF_NAME}_filtered_miss0_85.vcf.gz"

vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz" \
  --max-missing 0.80 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_80.vcf.gz"
cp "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_80.vcf.gz" \
   "${WORKDIR}/${FILT_DIR}/${REF_NAME}_filtered_miss0_80.vcf.gz"

vcftools --gzvcf "${WORKDIR}/${FILT_DIR}/filter_steps/4_maf.vcf.gz" \
  --max-missing 0.75 --recode --stdout \
  | gzip -c > "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_75.vcf.gz"
cp "${WORKDIR}/${FILT_DIR}/filter_steps/5_missing_0_75.vcf.gz" \
   "${WORKDIR}/${FILT_DIR}/${REF_NAME}_filtered_miss0_75.vcf.gz"
