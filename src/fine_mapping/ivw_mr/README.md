# Inverse Variance Weighted (IVW) Mendelian Randomization

- Run IVW TwosampleMR for cerebral cortical eQTLs against a given GWAS outcomes. 

## Step 1. Chunkize the input file into gene-level chunks

SUMSTATS: eQTL association data in summary statistics format
CHK_PREFIX: Prefix for chunk files
CHK_OUTDIR: Output directory for chunks

```
python3 chunkize_sumstats.py $SUMSTATS $CHK_PREFIX -o $CHK_OUTDIR
```

## Step 2. 

```
BFILE=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/genome/plink/snp_code_major/CerebralCortex.MetaAnalysis
EQTL=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/eqtl_discovery/sumstats/cortex_region_eQTLS.5e-8.joint.sumstats.txt
OUTCOME_FILE=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/outcome_data.txt
OUTDIR=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/ivw_mr
cd /gscmnt/gc2645/wgs/qtl/BrainQTL/src/new_src/ivw_mr
python3 ivw_mr_parallel.py $BFILE $EQTL  $OUTCOME_FILE -o $OUTDIR --chunksize 2
```
