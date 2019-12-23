# Formatting GWAS summary statistics for coloc and MR analysis

- The source GWAS summary statistics were formatted with the following commands
- This takes a large ammount of RAM so you need to request at least 75G on an interactive MGI computing blade 

```
OUTDIR=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data
KUNKLE=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/original/Kunkle_etal_Stage1_results.txt
RHEENEN=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/original/als.sumstats.lmm.All.txt
NALLS=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/original/METAL_PD_SE1_rs_sort_avsnp150_sorted_multianno.hg19_multianno_multianno_sorted_annovar_EAF_final.txt
cd /gscmnt/gc2645/wgs/qtl/BrainQTL/src/new_src/format

python3 format_gwas.py $KUNKLE kunkle AD_kunkle -o $OUTDIR
python3 format_gwas.py $RHEENEN rheenen ALS_rheenen -o $OUTDIR 
python3 format_gwas.py $NALLS metal_pd PD_nalls -o $OUTDIR 
```
