<p align="left">
<img src="assets/logo.png" alt="logo" width="200"/>
</p>

## Overview

**omniQTL** is a Python package developed at the [Translational Genomics Lab](https://med.stanford.edu/genomics-of-diabetes.html), led by Dr. Anna Gloyn at Stanford University. It is designed to enable robust, transparent, and reproducible quantitative trait loci (QTL) mapping across multiple molecular phenotypes, including caQTLs, eQTLs, pQTLs, and more.

The package provides a unified framework that spans the entire analysis pipeline — from quality control of genotyping and sequencing data, through genotype imputation, peak calling, gene quantification, and normalization, to QTL mapping (nominal, permute, and conditional) using QTLtools, followed by downstream analyses such as colocalization with (harmonized) GWAS summary statistics and visualization (e.g., locus zoom plots and genome tracks).

By integrating multi-level molecular phenotypes, omniQTL facilitates a deeper understanding of the genetic architecture underlying complex traits and diseases.

## Installation

```
git clone git@github.com:HaniceSun/omniQTL.git
conda env create -f environment.yaml
conda activate omniQTL
```

## Usage

```python
- genotyping

from omniQTL.genotyping import Genotyping
geno = Genotyping(bfile='Stanford')
geno.check_missingness()
geno.check_sex()
geno.check_heterozygosity()
geno.check_relatedness()
geno.check_hwe()
geno.pre_imputation_check()
geno.submit_to_topmed()
geno.get_topmed_results()

vcf_files = ['Stanford_TOPMed/Stanford_imputed.vcf.gz', 'Oxford_TOPMed/Oxford_imputed.vcf.gz']
output_vcf = 'Genotyping_imputed_OxfordStanford.vcf.gz'
geno.merge_vcfs(vcf_files, output_vcf)

geno.subset_samples_vcf(vcf_file='Genotyping_imputed_OxfordStanford.vcf.gz',
sample_list='caQTL_subset_sample_list.txt', output_vcf='caQTL_genotyping.vcf.gz')

geno.rename_samples_vcf(vcf_file='caQTL_genotyping.vcf.gz',
sample_mapping_file='caQTL_subset_sample_mapping.txt')

geno.annotate_variants_vcf(vcf_file='caQTL_genotyping_sampleRenamed.vcf.gz',
dbsnp_vcf='dbSNP157_GRCh38.vcf.gz', ref_fa='hg38.fa', norm_variants=False)

geno.filter_variants_vcf('caQTL_genotyping_sampleRenamed_rsID.vcf.gz')

- caQTL

from omniQTL.caQTL import CAQTL
omni = CAQTL()
omni.run_atac_seq_pipeline_encode(fq_file_list='bulkATACseqStanford.txt')
omni.get_pipeline_output()
omni.bam_flagstat(bam_dir='bams')
omni.get_numer_reads()
omni.get_percent_mapped_reads()
omni.get_tss_score()
omni.get_mbv_script()

omni.get_concensus_peaks(in_dir='peaks_qvalue', out_file='ATACseq_qvalue_consensus_peaks.bed')

omni.get_summit_extended_fixed_width_peaks(in_dir='peaks_qvalue',
out_file='ATACseq_qvalue_summitExtended_peaks.bed')

omni.get_peak_counts(in_file='ATACseq_qvalue_consensus_peaks.bed',
out_file='featureCounts_qvalue.sh')

omni.merge_counts_tables(counts_tables='counts_tables_qvalue.txt',
out_file='ATACseq_qvalue_peakCounts.txt')

omni.get_closest_genes(in_file='ATACseq_qvalue_peakCounts.txt',
gtf_file='Ensembl/Homo_sapiens.GRCh38.110.gtf')

omni.counts_to_tpm(in_file='ATACseq_qvalue_peakCounts_closestGene.txt')

omni.filter_phenotype_features(tpm_file='ATACseq_qvalue_peakCounts_closestGene_TPM.txt',
counts_file='ATACseq_qvalue_peakCounts_closestGene.txt')

omni.make_bed_for_QTLtools('ATACseq_qvalue_peakCounts_closestGene_TPM_peakFiltered.txt')
omni.run_PCA_on_bed(in_file='ATACseq_qvalue_peakCounts_closestGene_TPM_peakFiltered.bed.gz')
omni.get_QTLtools_script(pheno_file='ATACseq_qvalue_peakCounts_closestGene_TPM_peakFiltered.bed.gz', 

geno_file='caQTL_genotyping_sampleRenamed_rsID_variantFiltered.vcf.gz',
cov_file='ATACseq_qvalue_peakCounts_closestGene_TPM_peakFiltered_PC25.txt', out_suffix='qvalue')

- colocalization
from omniQTL.coloc import Coloc
cc = Coloc()
cc.download_gwas_harmoniser_reference()

```

## Author and License

**Author:** Han Sun

**Email:** hansun@stanford.edu

**License:** [MIT License](LICENSE)
